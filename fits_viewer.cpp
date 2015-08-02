#include "fits_viewer.h"

#include"fitsio.h"
#include"gsl/gsl_multifit.h"
#include"gsl/gsl_sort_vector.h"

#include<cstring>
#include<algorithm>
#include<random>

#include<QFile>
#include<QImage>
#include<QTransform>

#include<iostream>
#include<QDebug>


Fits_viewer::Fits_viewer(QString fits_filename, QWidget *parent): QGraphicsView(parent),
    currentPixmap(QPixmap()), currentCT(QVector<QRgb>(FITS_VIEWER_CT_LENGTH)), currentCursor(QCursor()),
    rubberBand_origin(QPoint()), rubberBand_end(QPoint()), rubberBand_pen(QPen("red"))
{
    gsl_set_error_handler_off();

    currentError = FITS_VIEWER_ERROR_OK;

    isImageLoaded = false;
    currentFITS_filename = fits_filename;
    currentImage_buffer = nullptr;
    currentImage_npix = 0;
    currentScaledImage_buffer = nullptr;

    _numScheduledScalings = 0;
    currentScale = 1.0;

    QObject::connect(this,SIGNAL(Fits_viewer_error(int)),this,SLOT(SetError(int)));

    scene = new QGraphicsScene(this);
    setScene(scene);

    GenerateCT(Fits_viewer::CT_NEGBW);
//    GenerateCT(Fits_viewer::CT_BW);

    QTransform tr(1.0,0.0,0.0,-1.0,0.0,0.0); // reflection about x-axis to put
    this->setTransform(tr);                  // the origin to bottom-left conner

    currentCursor.setShape(Qt::CrossCursor);
    this->setCursor(currentCursor);

    rubberBand = nullptr;
    isRubberBandActive = false;
    isRubberBandShown = false;


    rubberBand_pen.setWidth(0);

    pixmap_item = nullptr;

//    this->setDragMode(QGraphicsView::RubberBandDrag);
    this->setHorizontalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    this->setVerticalScrollBarPolicy(Qt::ScrollBarAsNeeded);
    this->setSizeAdjustPolicy(QAbstractScrollArea::AdjustToContentsOnFirstShow);
//    this->setInteractive(false);
//    this->setMouseTracking(false);

    this->setMouseTracking(true);

    currentFITS_filename = currentFITS_filename.trimmed();
    if ( currentFITS_filename.isEmpty() ) return;
    LoadFile(currentFITS_filename);
    isImageLoaded = true;

//    this->setMouseTracking(true);
//    this->setInteractive(true);

    double ch, cl;
    Compute_ZScale(&cl,&ch);
//    std::cout << "CUTS: " << cl << ", " << ch << "\n";

    ScaleImage(cl,ch);
}

Fits_viewer::Fits_viewer(QWidget *parent): Fits_viewer("", parent)
{
}


Fits_viewer::~Fits_viewer()
{
    delete[] currentImage_buffer;
    delete[] currentScaledImage_buffer;
}


        /*    Public methods    */

int Fits_viewer::getCurrentError() const
{
    return currentError;
}


void Fits_viewer::GetCurrentCuts(double *low_cuts, double *high_cuts)
{
    *low_cuts = currentLowCut;
    *high_cuts = currentHighCut;
}


void Fits_viewer::GetImageMinMax(double *min_val, double *max_val)
{
    *min_val = currentImageMinVal;
    *max_val = currentImageMaxVal;
}

        /*    Public slots    */

void Fits_viewer::LoadFile(QString fits_filename)
{
    QString str = fits_filename.trimmed();
    if ( str.isEmpty() || str.isNull() ) return;

    fitsfile *FITS_fptr;

    int fits_status = 0;
    currentError = FITS_VIEWER_ERROR_OK;

    // THE ONLY 2D-images is support now!!!
    int maxdim = 2;
    long naxes[maxdim];
    int naxis, bitpix;
    LONGLONG nelem = 1;


    char* filename = fits_filename.toLocal8Bit().data();

    try {
        fits_open_image(&FITS_fptr, filename, READONLY, &fits_status);
        if ( fits_status ) throw fits_status;

        fits_read_imghdr(FITS_fptr, maxdim, NULL, &bitpix, &naxis, naxes, NULL, NULL, NULL, &fits_status);
        if ( fits_status ) throw fits_status;

        for ( int i = 0; i < maxdim; ++i ) {
            nelem *= naxes[i];
            currentImage_dim[i] = naxes[i];
        }

        currentImage_npix = nelem;
        delete[] currentImage_buffer;
        currentImage_buffer = new double[currentImage_npix];

        fits_read_img(FITS_fptr, TDOUBLE, 1, nelem, NULL, (void*) currentImage_buffer, NULL, &fits_status);
        if ( fits_status ) throw fits_status;

        fits_close_file(FITS_fptr, &fits_status);
        if ( fits_status ) throw fits_status;

    } catch (std::bad_alloc &ex) {
        currentImage_buffer = nullptr;
        throw (int)FITS_VIEWER_ERROR_BAD_ALLOC;
    } catch (int err) {
        emit Fits_viewer_error(err);
        fits_close_file(FITS_fptr, &fits_status);

        delete[] currentImage_buffer;
        currentImage_buffer = nullptr;
        currentImage_npix = 0;

        return;
    }

    auto minmax = std::minmax_element(currentImage_buffer,currentImage_buffer+currentImage_npix);
    currentImageMinVal = *minmax.first;
    currentImageMaxVal = *minmax.second;

//    this->setInteractive(true);
//    this->setMouseTracking(true);

    currentScale = 1.0;

    Compute_ZScale(&currentLowCut,&currentHighCut);
    ScaleImage(currentLowCut,currentHighCut);
//    ShowImage();
    isImageLoaded = true;
}


void Fits_viewer::ScaleImage(const double low_val, const double high_val)
{
    double scaled_val;
    double low_cuts = low_val;
    double high_cuts = high_val;

    if ( (currentImage_buffer == nullptr) || (currentImage_npix == 0) ) return;

    if ( low_val >= high_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    if ( currentImageMinVal > high_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    if ( currentImageMaxVal < low_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    if ( low_val < currentImageMinVal ) low_cuts = currentImageMinVal;
    if ( high_val > currentImageMaxVal ) high_cuts = currentImageMaxVal;

    double range = high_cuts-low_cuts;

    try {
        delete[] currentScaledImage_buffer;
        currentScaledImage_buffer = new uchar[currentImage_npix];

        for ( size_t i = 0; i < currentImage_npix; ++i ) {
            if ( currentImage_buffer[i] <= low_val ) {
                currentScaledImage_buffer[i] = 0;
                continue;
            }
            if ( currentImage_buffer[i] >= high_val ) {
                currentScaledImage_buffer[i] = FITS_VIEWER_CT_LENGTH-1;
                continue;
            }
            scaled_val = (currentImage_buffer[i]-low_val)/range;
            currentScaledImage_buffer[i] = static_cast<uchar>(std::lround(scaled_val*255));
        }

        currentLowCut = low_cuts;
        currentHighCut = high_cuts;

        emit ScalingIsChanged(currentLowCut,currentHighCut);

        ShowImage();
    } catch (std::bad_alloc &ex) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_ALLOC);
        return;
    }

}


/*    Protected slots   */

void Fits_viewer::mouseMoveEvent(QMouseEvent* event)
{
    if ( !isImageLoaded ) return; // no image

    QPointF pos =  mapToScene( event->pos() );

    if ( (pos.x() >= 0) && (pos.x() < currentImage_dim[0]) &&
         (pos.y() >= 0) && (pos.y() < currentImage_dim[1]) ) {
        uint x,y;
        x = (uint)pos.x();
        y = (uint)pos.y();

        double val = currentImage_buffer[y*currentImage_dim[0]+x];

        // convert to FITS pixel coordinate notation (starting from 1, integer coordinates are at the pixel's center)
        pos.setX(pos.x()+0.5);
        pos.setY(pos.y()+0.5);

        emit ImagePoint(pos,val);
        emit ImagePoint(pos.x(),pos.y(),val);
    }

    if ( isRubberBandActive && (event->buttons() & Qt::LeftButton) ) {

        rubberBand->setVisible(true);
        isRubberBandShown = true;

        rubberBand_end = mapToScene(event->pos());

        if ( rubberBand_end.x() < 0 ) {
            rubberBand_end.setX(0.0);
        }
        if ( rubberBand_end.y() < 0 ) {
            rubberBand_end.setY(0.0);
        }
        if ( rubberBand_end.x() >= currentImage_dim[0] ) {
            rubberBand_end.setX(currentImage_dim[0]-1);
        }
        if ( rubberBand_end.y() >= currentImage_dim[1] ) {
            rubberBand_end.setY(currentImage_dim[1]-1);
        }

        rubberBand->setRect(QRectF(rubberBand_origin, rubberBand_end).normalized());
    }
}


void Fits_viewer::mousePressEvent(QMouseEvent *event)
{
    if ( event->button() == Qt::LeftButton ) {
        if ( isRubberBandShown ) {
            scene->removeItem(rubberBand);
            isRubberBandShown = false;
        }

        rubberBand_origin = mapToScene(event->pos());
//        qDebug() << rubberBand_origin;

        // prevent rectangle conner is out of image
        if ( rubberBand_origin.x() < 0 ) {
            rubberBand_origin.setX(0.0);
        }
        if ( rubberBand_origin.x() >= currentImage_dim[0] ) {
            rubberBand_origin.setX(currentImage_dim[0]-1);
        }

        if ( rubberBand_origin.y() < 0 ) {
            rubberBand_origin.setY(0.0);
        }
        if ( rubberBand_origin.y() >= currentImage_dim[1] ) {
            rubberBand_origin.setY(currentImage_dim[1]-1);
        }

        rubberBand_end = QPointF(rubberBand_origin);
        rubberBand = scene->addRect(QRectF(rubberBand_origin, QSize()),rubberBand_pen);

        rubberBand->setVisible(false);
        isRubberBandActive = true;
    }

    if ( event->button() == Qt::RightButton ) {
        if ( isRubberBandShown ) {
            scene->removeItem(rubberBand);
            isRubberBandActive = false;
            isRubberBandShown = false;

            double zmin, zmax;

            Compute_ZScale(&zmin,&zmax);
            ScaleImage(zmin,zmax);
//            ShowImage();
//            std::cout << "CUTS: [" << zmin << ", " << zmax << "]\n";
        }
    }
}


void Fits_viewer::mouseReleaseEvent(QMouseEvent *event)
{
    if ( event->button() == Qt::LeftButton ) {
        if ( isRubberBandShown ) {
            isRubberBandActive = false;
            QRectF rect = rubberBand->rect();
//            std::cout << "SELECT: " << rect.x() << ", " << rect.y() << std::endl;
            emit SelectedRegion(rect);
        }
    }
}


void Fits_viewer::resizeEvent(QResizeEvent *event)
{
    if ( pixmap_item != nullptr ) {
        this->fitInView(pixmap_item,Qt::KeepAspectRatio);
        this->scale(currentScale, currentScale);
    }
}


void Fits_viewer::wheelEvent(QWheelEvent *event)
{
    int numDegrees = event->delta() / 8;
    int numSteps = numDegrees / 15; // see QWheelEvent documentation
    _numScheduledScalings += numSteps;
    if (_numScheduledScalings * numSteps < 0) // if user moved the wheel in another direction, we reset previously scheduled scalings
        _numScheduledScalings = numSteps;

    QTimeLine *anim = new QTimeLine(350, this);
    anim->setUpdateInterval(20);

    connect(anim, SIGNAL (valueChanged(qreal)), SLOT (scalingTime(qreal)));
    connect(anim, SIGNAL (finished()), SLOT (animFinished()));
    anim->start();
}


void Fits_viewer::keyPressEvent(QKeyEvent *event)
{
//    std::cout << "key pressed: " << std::hex << event->key() << "\n";
    qreal factor;
    bool scaling = false;
    if ( event->key() == Qt::Key_1 ) {
        factor = 1.0;
        scaling = true;
    }
    if ( event->key() == Qt::Key_2 ) {
        factor = 2.0;
        scaling = true;
    }
    if ( event->key() == Qt::Key_3 ) {
        factor = 3.0;
        scaling = true;
    }
    if ( event->key() == Qt::Key_4 ) {
        factor = 4.0;
        scaling = true;
    }
    if ( event->key() == Qt::Key_Escape ) {
        factor = 1.0;
        scaling = true;
    }
    if ( scaling ) {
        qreal tmp = factor;
        factor /= currentScale;
        currentScale = tmp;
        this->scale(factor,factor);
    }
}


/*    Private slots   */

void Fits_viewer::scalingTime(qreal x)
{
    qreal factor = 1.0+ qreal(_numScheduledScalings) / 300.0;
    currentScale *= factor;
    if ( currentScale <= 1.0 ) {
        currentScale = 1.0;
        return;
    }
    scale(factor, factor);
}


void Fits_viewer::animFinished()
{
    if (_numScheduledScalings > 0) _numScheduledScalings--; else _numScheduledScalings++;
    sender()->~QObject();
}

void Fits_viewer::SetError(int err)
{
    currentError = err;
}


/*    Private methods   */

void Fits_viewer::ShowImage()
{

    QImage im = QImage(currentScaledImage_buffer,currentImage_dim[0],currentImage_dim[1],currentImage_dim[0],QImage::Format_Indexed8);
    im.setColorTable(currentCT);

    currentPixmap = QPixmap::fromImage(im);

    scene->clear();
    scene->setSceneRect(0,0,currentImage_dim[0],currentImage_dim[1]);
//    qDebug() << scene->sceneRect();

//    pixmap_item = scene->addPixmap(currentPixmap.scaledToWidth(this->width()));
    pixmap_item = scene->addPixmap(currentPixmap);
    this->centerOn(pixmap_item);
    this->fitInView(pixmap_item,Qt::KeepAspectRatio);
    this->scale(currentScale,currentScale);
//    pixmap_item->setFlag(QGraphicsItem::ItemIsMovable);
}


void Fits_viewer::GenerateCT(Fits_viewer::ColorTable ct)
{
    int j;
    qreal ct_step = 255.0/(FITS_VIEWER_CT_LENGTH-1);

    switch (ct) {
        case Fits_viewer::CT_BW: { // black-and-white (grayscale)
            for ( int i = 0; i < FITS_VIEWER_CT_LENGTH; ++i ) {
                j = i*ct_step;
                currentCT[i] = qRgb(j,j,j);
//                std::cout << "color = " << j << std::endl;
            }
            break;
        }
        case Fits_viewer::CT_NEGBW: { // negative grayscale
            for ( int i = 0; i < FITS_VIEWER_CT_LENGTH; ++i ) {
                j = 255-i*ct_step;
                if ( j < 0 ) j = 0;
//                std::cout << "color = " << j << std::endl;
                currentCT[i] = qRgb(j,j,j);
            }
            break;
        }
        default: emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_COLOR_TABLE);
    }
}


void Fits_viewer::Compute_ZScale(double *zmin, double *zmax, const double contrast)
{
    // contrast parameter must be in range of (0,1]
    double cst = (contrast <= 0.0) ? 0.25 : contrast;
    cst = (cst > 1.0) ? 1.0 : cst;

    bool max_nsample = false;
    size_t nsample, xl, xr, yl, yr;
    size_t midpoint;

    gsl_vector *yy, *cc;
    gsl_matrix *cov, *X;
    gsl_multifit_robust_workspace *work;

    QRectF rect = QRectF(rubberBand_origin, rubberBand_end).normalized();

    xl = std::roundl(rect.topLeft().x());
    xr = std::roundl(rect.bottomRight().x());
    yl = std::roundl(rect.topLeft().y());
    yr = std::roundl(rect.bottomRight().y());

//    std::cout << "xl = " << xl << ", xr = " << xr << ", yl = " << yl << ", yr = " << yr << "\n";

    if ( xl == xr ) {
        xl = 0;
        xr = currentImage_dim[0]-1;
    }
    if ( yl == yr ) {
        yl = 0;
        yr = currentImage_dim[1]-1;
    }

    nsample = (xr-xl+1)*(yr-yl+1);

    if ( nsample > FITS_VIEWER_ZSCALE_NSAMPLE ) {
        nsample = FITS_VIEWER_ZSCALE_NSAMPLE;
        max_nsample = true;
    }


//    std::cout << "sub-im: [" << xl << ":" << xr << "," << yl << ":" << yr << "]\n";
//    std::cout << "nsample: " << nsample << std::endl;

    yy = NULL;
    X = NULL;
    work = NULL;
    cov = NULL;
    cc = NULL;

    midpoint = nsample/2;

    try {
        // copy sub-image

        yy = gsl_vector_alloc(nsample);
        if ( yy == NULL ) throw 10;

        X = gsl_matrix_alloc(nsample,2);
        if ( X == NULL ) throw 10;

        work = gsl_multifit_robust_alloc (gsl_multifit_robust_bisquare, X->size1, X->size2);
        if ( work == NULL ) throw 10;

        cc = gsl_vector_alloc(2);
        if ( cc == NULL ) throw 10;

        cov = gsl_matrix_alloc(2,2);
        if ( cov == NULL ) throw 10;

        size_t j = 0;

        if ( max_nsample ) { // generate random sample coordinates
            size_t x,y;
            std::random_device rd;
            std::mt19937 gen(rd());
            std::uniform_int_distribution<size_t> disX(xl, xr);
            std::uniform_int_distribution<size_t> disY(yl, yr);
            for ( size_t i = 0; i < nsample; ++i ) {
                x = disX(gen);
                y = disY(gen);
                gsl_vector_set(yy,i,currentImage_buffer[y*currentImage_dim[0]+x]);
            }
        } else {
            for ( size_t y = yl; y <= yr; ++y ) {
                for ( size_t x = xl; x <= xr; ++x ) {
                    gsl_vector_set(yy,j,currentImage_buffer[y*currentImage_dim[0]+x]);
                    ++j;
                }
            }
        }

        gsl_sort_vector(yy);

        // fill matrix
        for ( size_t i = 0; i < nsample; ++i ) {
            gsl_matrix_set(X,i,0,1.0);
            gsl_matrix_set(X,i,1,(double)i);
        }

        double zmed;
        if ( (nsample % 2) == 1 ) {
            zmed = gsl_vector_get(yy,midpoint);
        } else {
            zmed = (gsl_vector_get(yy,midpoint) + gsl_vector_get(yy,midpoint-1))/2.0;
        }

        double data_min = gsl_vector_get(yy,0);
        double data_max = gsl_vector_get(yy,nsample-1);

        // robust linear fit
        int code = gsl_multifit_robust (X, yy, cc, cov, work);

        if ( !code ) {
//            std::cout << "midpoint: " << midpoint << "\n";
//            std::cout << "zmed: " << zmed << ", contrast: " << cst << "\n";
//            std::cout << "coeffs: " << gsl_vector_get(cc,0) << ", " << gsl_vector_get(cc,1) << "\n";
            *zmin = zmed + (gsl_vector_get(cc,1)/cst) * (1.0 - midpoint);
            *zmax = zmed + (gsl_vector_get(cc,1)/cst) * (1.0*nsample - midpoint);

            if ( *zmin < data_min ) *zmin = data_min;
            if ( *zmax > data_max ) *zmax = data_max;
        } else {
            *zmin = data_min;
            *zmax = data_max;
        }

    } catch (int err) {
        gsl_vector_free(yy);
        gsl_matrix_free(X);
        gsl_multifit_robust_free(work);
        gsl_matrix_free(cov);
        gsl_vector_free(cc);

        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_ALLOC);
        return;
    }

    gsl_vector_free(yy);
    gsl_matrix_free(X);
    gsl_multifit_robust_free(work);
    gsl_matrix_free(cov);
    gsl_vector_free(cc);
}
