TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle
CONFIG -= qt

SOURCES += main.cpp

HEADERS += \
    dynamicfloat.h

QMAKE_CXXFLAGS += -std=c++14