#include "fits_viewer.h"

#include"fitsio.h"

#include<algorithm>

#include<QFile>
#include<QImage>


#include<iostream>


Fits_viewer::Fits_viewer(QString fits_filename, QWidget *parent): QGraphicsView(parent),
    currentPixmap(QPixmap()), currentCT(QVector<QRgb>(256)), currentCursor(QCursor())
{
    currentError = 0;

    currentFITS_filename = fits_filename;
    currentImage_buffer = nullptr;
    currentImage_npix = 0;
    currentScaledImage_buffer = nullptr;

    QObject::connect(this,SIGNAL(Fits_viewer_error(int)),this,SLOT(SetError(int)));

    scene = new QGraphicsScene(this);
    setScene(scene);

    GenerateCT(Fits_viewer::CT_NEGBW);

    LoadFile(currentFITS_filename);
    ScaleImage(100,1000);
    ShowImage();

    currentCursor.setShape(Qt::CrossCursor);
    this->setCursor(currentCursor);

    this->setMouseTracking(true);
    this->setDragMode(QGraphicsView::NoDrag);
    this->setHorizontalScrollBarPolicy(Qt::ScrollBarAlwaysOff);
    this->setVerticalScrollBarPolicy(Qt::ScrollBarAlwaysOff);

    currentScale = 1.0;
}

Fits_viewer::Fits_viewer(QWidget *parent): Fits_viewer("", parent)
{
}


Fits_viewer::~Fits_viewer()
{
    delete[] currentImage_buffer;
    delete[] currentScaledImage_buffer;
}


/*    Public slots    */

void Fits_viewer::LoadFile(QString fits_filename)
{
    fitsfile *FITS_fptr;

    int fits_status = 0;

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

}


void Fits_viewer::ScaleImage(const double low_val, const double high_val)
{
    double scaled_val;

    if ( (currentImage_buffer == nullptr) || (currentImage_npix == 0) ) return;

    if ( low_val >= high_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    auto minmax = std::minmax_element(currentImage_buffer,currentImage_buffer+currentImage_npix);

    if ( *minmax.first > high_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    if ( *minmax.second < low_val ) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_SCALE_PARS);
        return;
    }

    try {
        delete[] currentScaledImage_buffer;
        currentScaledImage_buffer = new uchar[currentImage_npix];

        for ( size_t i = 0; i < currentImage_npix; ++i ) {
            if ( currentImage_buffer[i] <= low_val ) {
                currentScaledImage_buffer[i] = 0;
                continue;
            }
            if ( currentImage_buffer[i] >= high_val ) {
                currentScaledImage_buffer[i] = 255;
                continue;
            }
            scaled_val = static_cast<uchar>(std::lround((currentImage_buffer[i]-low_val)/(*minmax.second-*minmax.first)));
            currentScaledImage_buffer[i] = scaled_val*255;
        }
    } catch (std::bad_alloc &ex) {
        emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_ALLOC);
        return;
    }

}


/*    Protected slots   */

void Fits_viewer::mouseMoveEvent(QMouseEvent* event)
{
    QPointF pos =  mapToScene( event->pos() );
//    std::cout << "POS: [" << pos.x() << ": " << pos.y() << "]  ";
    if ( (pos.x() >= 0) && (pos.x() < currentImage_dim[0]) &&
         (pos.y() >= 0) && (pos.y() < currentImage_dim[1]) ) {
        uint x,y;
        x = (uint)pos.x();
        y = (uint)pos.y();
//        std::cout << "VAL: " << currentImage_buffer[y*currentImage_dim[0]+x] << std::endl;
    }
}


void Fits_viewer::resizeEvent(QResizeEvent *event)
{
    this->fitInView(pixmap_item,Qt::KeepAspectRatio);
    this->scale(currentScale, currentScale);
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




/*    Private slots   */

void Fits_viewer::scalingTime(qreal x)
{
    qreal factor = 1.0+ qreal(_numScheduledScalings) / 300.0;
    scale(factor, factor);
    currentScale *= factor;

    std::cout << "scale: " << factor << std::endl;
    std::cout << "scling: " << _numScheduledScalings << std::endl;
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
//    pixmap_item = scene->addPixmap(currentPixmap.scaledToWidth(this->width()));
    pixmap_item = scene->addPixmap(currentPixmap);
    this->fitInView(pixmap_item,Qt::KeepAspectRatio);
//    pixmap_item->setFlag(QGraphicsItem::ItemIsMovable);
}


void Fits_viewer::GenerateCT(Fits_viewer::ColorTable ct)
{
    switch (ct) {
        case Fits_viewer::CT_BW: { // black-and-white (grayscale)
            for ( int i = 0; i < 256; ++i ) currentCT[i] = qRgb(i,i,i);
        }
        case Fits_viewer::CT_NEGBW: { // negative grayscale
            int j;
            for ( int i = 0; i < 256; ++i ) {
                j = 255-i;
                currentCT[i] = qRgb(j,j,j);
            }
        }
        default: emit Fits_viewer_error(FITS_VIEWER_ERROR_BAD_COLOR_TABLE);
    }
}
