#-------------------------------------------------
#
# Project created by QtCreator 2015-06-06T17:07:01
#
#-------------------------------------------------

QT       += widgets

TARGET = fits_viewer
TEMPLATE = lib

DEFINES += FITS_VIEWER_LIBRARY

SOURCES += fits_viewer.cpp

HEADERS += fits_viewer.h\
        fits_viewer_global.h

unix {
    target.path = /usr/lib
    INSTALLS += target
}
