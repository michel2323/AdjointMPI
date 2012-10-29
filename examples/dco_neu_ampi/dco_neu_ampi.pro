#-------------------------------------------------
#
# Project created by QtCreator 2011-05-02T12:41:53
#
#-------------------------------------------------

QT       -= core gui

TARGET = dco_neu_ampi
CONFIG   += console
CONFIG   -= app_bundle

TEMPLATE = app


SOURCES += main.cpp

SOURCES += dco_tape.cpp
HEADERS += dco_tape.cpp

QMAKE_CXX = mpic++

INCLUDEPATH += ../../include


LIBS += ../../src/libAMPI.a -lmpi
