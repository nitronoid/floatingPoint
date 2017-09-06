
QT += testlib
QT -= gui

TARGET = main

TEMPLATE = app
CONFIG += console c++14
CONFIG -= app_bundle

SOURCES += \
    main.cpp \
    dynamicfloat.inl

HEADERS += \
    dynamicfloat.h

QMAKE_CXXFLAGS += -std=c++14

#clang
linux-clang++: QMAKE_CXXFLAGS += -Weverything -Wno-c++98-compat
#gcc
linux-g++: QMAKE_CXXFLAGS += -Wall -Wextra -pedantic-errors

DEFINES += SRCDIR=\\\"$$PWD/\\\"
