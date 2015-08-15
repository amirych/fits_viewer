#ifndef FITS_VIEWER_H
#define FITS_VIEWER_H

#include "fits_viewer_global.h"

#include<QtWidgets>
#include<QPixmap>
#include<QVector>
#include<QRgb>
#include<QCursor>
#include<QRubberBand>
#include<QGraphicsRectItem>

#define FITS_VIEWER_ERROR_OK 0
#define FITS_VIEWER_ERROR_BAD_ALLOC 10000
#define FITS_VIEWER_ERROR_BAD_SCALE_PARS 10001
#define FITS_VIEWER_ERROR_BAD_COLOR_TABLE 10002

#define FITS_VIEWER_CT_LENGTH 256

#define FITS_VIEWER_ZSCALE_NSAMPLE 10000

class FITS_VIEWERSHARED_EXPORT Fits_viewer: public QGraphicsView
{

Q_OBJECT

public:
    enum ColorTable {CT_BW, CT_NEGBW};

    Fits_viewer(QString fits_filename, QWidget *parent = 0);
    Fits_viewer(QWidget *parent = 0);
    ~Fits_viewer();

    void GetCurrentCuts(double *low_cuts, double *high_cuts);
    void GetImageMinMax(double *min_val, double *max_val);
    int getCurrentError() const;

    double* getCurrentSubImage(QRect *rect);

signals:
    void Fits_viewer_error(int err);
    void ImagePoint(QPointF pos, double val);
    void ImagePoint(double xpos, double ypos, double val);
    void SelectedRegion(QRectF region);
    void DeselectRegion();
    void ScalingIsChanged(double low_cuts, double high_cuts);

public slots:
    void LoadFile(QString fits_filename);
    void ScaleImage(const double low_val, const double high_val);

protected:
    virtual void mouseMoveEvent(QMouseEvent* event);
    virtual void resizeEvent(QResizeEvent* event);
    virtual void wheelEvent(QWheelEvent* event);
    virtual void mousePressEvent(QMouseEvent* event);
    virtual void mouseReleaseEvent(QMouseEvent* event);
    virtual void keyPressEvent(QKeyEvent* event);


private slots:
    void SetError(int err);
    void scalingTime(qreal x);
    void animFinished();

private:
    int currentError;
    QString currentFITS_filename;

    bool isImageLoaded;
    double *currentImage_buffer;
    uchar *currentScaledImage_buffer;
    size_t currentImage_npix;
    size_t currentImage_dim[2];
    double currentImageMinVal;
    double currentImageMaxVal;
    double currentLowCut,currentHighCut;

    double *currentSubImage;


    QPixmap currentPixmap;
    QVector<QRgb> currentCT;
    QGraphicsScene *scene;
    QGraphicsItem *pixmap_item;
    QCursor currentCursor;
//    QRubberBand *rubberBand;
    QGraphicsRectItem *rubberBand;

    void ShowImage();
    void GenerateCT(Fits_viewer::ColorTable ct);
    void Compute_ZScale(double *zmin, double *zmax, const double contrast = 0.25);
    void ZScale(double *zmin, double *zmax, const double contrast = 0.25);

    int _numScheduledScalings;
    qreal currentScale;

    QPointF rubberBand_origin, rubberBand_end;
    QPen rubberBand_pen;
    bool isRubberBandActive;
    bool isRubberBandShown;
};

#endif // FITS_VIEWER_H
