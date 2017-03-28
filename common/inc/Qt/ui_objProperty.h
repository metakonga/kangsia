/********************************************************************************
** Form generated from reading UI file 'objPropertyH31192.ui'
**
** Created by: Qt User Interface Compiler version 5.4.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef OBJPROPERTYH31192_H
#define OBJPROPERTYH31192_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QDialog>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_ObjProperty
{
public:
    QGridLayout *gridLayout_2;
    QScrollArea *scrollArea;
    QWidget *scrollAreaWidgetContents;
    QGridLayout *gridLayout;
    QScrollArea *scrollArea_2;
    QWidget *scrollAreaWidgetContents_2;

    void setupUi(QDialog *ObjProperty)
    {
        if (ObjProperty->objectName().isEmpty())
            ObjProperty->setObjectName(QStringLiteral("ObjProperty"));
        ObjProperty->resize(400, 443);
        gridLayout_2 = new QGridLayout(ObjProperty);
        gridLayout_2->setObjectName(QStringLiteral("gridLayout_2"));
        scrollArea = new QScrollArea(ObjProperty);
        scrollArea->setObjectName(QStringLiteral("scrollArea"));
        scrollArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 380, 210));
        gridLayout = new QGridLayout(scrollAreaWidgetContents);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        scrollArea->setWidget(scrollAreaWidgetContents);

        gridLayout_2->addWidget(scrollArea, 0, 0, 1, 1);

        scrollArea_2 = new QScrollArea(ObjProperty);
        scrollArea_2->setObjectName(QStringLiteral("scrollArea_2"));
        scrollArea_2->setWidgetResizable(true);
        scrollAreaWidgetContents_2 = new QWidget();
        scrollAreaWidgetContents_2->setObjectName(QStringLiteral("scrollAreaWidgetContents_2"));
        scrollAreaWidgetContents_2->setGeometry(QRect(0, 0, 380, 205));
        scrollArea_2->setWidget(scrollAreaWidgetContents_2);

        gridLayout_2->addWidget(scrollArea_2, 1, 0, 1, 1);


        retranslateUi(ObjProperty);

        QMetaObject::connectSlotsByName(ObjProperty);
    } // setupUi

    void retranslateUi(QDialog *ObjProperty)
    {
        ObjProperty->setWindowTitle(QApplication::translate("ObjProperty", "Properies", 0));
    } // retranslateUi

};

namespace Ui {
    class ObjProperty: public Ui_ObjProperty {};
} // namespace Ui

QT_END_NAMESPACE

#endif // OBJPROPERTYH31192_H
