#ifndef FITS_VIEWER_H
#define FITS_VIEWER_H

#include "fits_viewer_global.h"

#include<QtWidgets>
#include<QPixmap>
#include<QVector>
#include<QRgb>
#include<QCursor>


#define FITS_VIEWER_ERROR_BAD_ALLOC 10000
#define FITS_VIEWER_ERROR_BAD_SCALE_PARS 10001
#define FITS_VIEWER_ERROR_BAD_COLOR_TABLE 10002

class FITS_VIEWERSHARED_EXPORT Fits_viewer: public QGraphicsView
{

Q_OBJECT

public:
    enum ColorTable {CT_BW, CT_NEGBW};

    Fits_viewer(QString fits_filename, QWidget *parent = 0);
    Fits_viewer(QWidget *parent = 0);
    ~Fits_viewer();

signals:
    void Fits_viewer_error(int err);

public slots:
    void LoadFile(QString fits_filename);
    void ScaleImage(const double low_val, const double high_val);

protected slots:
    void mouseMoveEvent(QMouseEvent* event);
    void resizeEvent(QResizeEvent* event);
    void wheelEvent(QWheelEvent* event);

private slots:
    void SetError(int err);
    void scalingTime(qreal x);
    void animFinished();

private:
    int currentError;
    QString currentFITS_filename;

    double *currentImage_buffer;
    uchar *currentScaledImage_buffer;
    size_t currentImage_npix;
    size_t currentImage_dim[2];

    QPixmap currentPixmap;
    QVector<QRgb> currentCT;
    QGraphicsScene *scene;
    QGraphicsItem *pixmap_item;
    QCursor currentCursor;

    void ShowImage();
    void GenerateCT(Fits_viewer::ColorTable ct);

    int _numScheduledScalings;
    qreal currentScale;
};

#endif // FITS_VIEWER_H
