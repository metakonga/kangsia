/********************************************************************************
** Form generated from reading UI file 'parview.ui'
**
** Created by: Qt User Interface Compiler version 5.4.2
**
** WARNING! All changes made in this file will be lost when recompiling UI file!
********************************************************************************/

#ifndef UI_PARVIEW_H
#define UI_PARVIEW_H

#include <QtCore/QVariant>
#include <QtWidgets/QAction>
#include <QtWidgets/QApplication>
#include <QtWidgets/QButtonGroup>
#include <QtWidgets/QGridLayout>
#include <QtWidgets/QHeaderView>
#include <QtWidgets/QMainWindow>
#include <QtWidgets/QMenu>
#include <QtWidgets/QMenuBar>
#include <QtWidgets/QScrollArea>
#include <QtWidgets/QStatusBar>
#include <QtWidgets/QToolBar>
#include <QtWidgets/QWidget>

QT_BEGIN_NAMESPACE

class Ui_parVIEW
{
public:
    QAction *actionChange_Shape;
    QAction *actionMilkShape_3D_ASCII;
    QAction *actionMBD_Result_ASCII;
    QAction *actionDEM_Result_ASCII;
    QAction *actionProperty;
    QWidget *centralWidget;
    QGridLayout *gridLayout;
    QScrollArea *GraphicArea;
    QWidget *scrollAreaWidgetContents;
    QMenuBar *menuBar;
    QMenu *menu;
    QMenu *menuImport;
    QMenu *menuExport;
    QToolBar *mainToolBar;
    QStatusBar *statusBar;
    QToolBar *secToolBar;

    void setupUi(QMainWindow *parVIEW)
    {
        if (parVIEW->objectName().isEmpty())
            parVIEW->setObjectName(QStringLiteral("parVIEW"));
        parVIEW->resize(600, 400);
        actionChange_Shape = new QAction(parVIEW);
        actionChange_Shape->setObjectName(QStringLiteral("actionChange_Shape"));
        actionChange_Shape->setCheckable(true);
        actionMilkShape_3D_ASCII = new QAction(parVIEW);
        actionMilkShape_3D_ASCII->setObjectName(QStringLiteral("actionMilkShape_3D_ASCII"));
        actionMBD_Result_ASCII = new QAction(parVIEW);
        actionMBD_Result_ASCII->setObjectName(QStringLiteral("actionMBD_Result_ASCII"));
        actionDEM_Result_ASCII = new QAction(parVIEW);
        actionDEM_Result_ASCII->setObjectName(QStringLiteral("actionDEM_Result_ASCII"));
        actionProperty = new QAction(parVIEW);
        actionProperty->setObjectName(QStringLiteral("actionProperty"));
        centralWidget = new QWidget(parVIEW);
        centralWidget->setObjectName(QStringLiteral("centralWidget"));
        gridLayout = new QGridLayout(centralWidget);
        gridLayout->setSpacing(6);
        gridLayout->setContentsMargins(11, 11, 11, 11);
        gridLayout->setObjectName(QStringLiteral("gridLayout"));
        GraphicArea = new QScrollArea(centralWidget);
        GraphicArea->setObjectName(QStringLiteral("GraphicArea"));
        GraphicArea->setWidgetResizable(true);
        scrollAreaWidgetContents = new QWidget();
        scrollAreaWidgetContents->setObjectName(QStringLiteral("scrollAreaWidgetContents"));
        scrollAreaWidgetContents->setGeometry(QRect(0, 0, 580, 315));
        GraphicArea->setWidget(scrollAreaWidgetContents);

        gridLayout->addWidget(GraphicArea, 0, 0, 1, 1);

        parVIEW->setCentralWidget(centralWidget);
        menuBar = new QMenuBar(parVIEW);
        menuBar->setObjectName(QStringLiteral("menuBar"));
        menuBar->setGeometry(QRect(0, 0, 600, 21));
        menu = new QMenu(menuBar);
        menu->setObjectName(QStringLiteral("menu"));
        menuImport = new QMenu(menu);
        menuImport->setObjectName(QStringLiteral("menuImport"));
        menuExport = new QMenu(menu);
        menuExport->setObjectName(QStringLiteral("menuExport"));
        parVIEW->setMenuBar(menuBar);
        mainToolBar = new QToolBar(parVIEW);
        mainToolBar->setObjectName(QStringLiteral("mainToolBar"));
        mainToolBar->setMovable(false);
        mainToolBar->setFloatable(false);
        parVIEW->addToolBar(Qt::TopToolBarArea, mainToolBar);
        statusBar = new QStatusBar(parVIEW);
        statusBar->setObjectName(QStringLiteral("statusBar"));
        parVIEW->setStatusBar(statusBar);
        secToolBar = new QToolBar(parVIEW);
        secToolBar->setObjectName(QStringLiteral("secToolBar"));
        parVIEW->addToolBar(Qt::TopToolBarArea, secToolBar);
        parVIEW->insertToolBarBreak(secToolBar);

        menuBar->addAction(menu->menuAction());
        menu->addAction(actionChange_Shape);
        menu->addAction(menuExport->menuAction());
        menu->addAction(menuImport->menuAction());
        menu->addAction(actionProperty);
        menuImport->addAction(actionMilkShape_3D_ASCII);
        menuExport->addAction(actionMBD_Result_ASCII);
        menuExport->addAction(actionDEM_Result_ASCII);

        retranslateUi(parVIEW);

        QMetaObject::connectSlotsByName(parVIEW);
    } // setupUi

    void retranslateUi(QMainWindow *parVIEW)
    {
        parVIEW->setWindowTitle(QApplication::translate("parVIEW", "parVIEW", 0));
        actionChange_Shape->setText(QApplication::translate("parVIEW", "Change Shape", 0));
        actionMilkShape_3D_ASCII->setText(QApplication::translate("parVIEW", "MilkShape 3D ASCII", 0));
        actionMBD_Result_ASCII->setText(QApplication::translate("parVIEW", "MBD Result ASCII", 0));
        actionDEM_Result_ASCII->setText(QApplication::translate("parVIEW", "DEM Result ASCII", 0));
        actionProperty->setText(QApplication::translate("parVIEW", "Property", 0));
        menu->setTitle(QApplication::translate("parVIEW", "\352\270\260\353\212\245", 0));
        menuImport->setTitle(QApplication::translate("parVIEW", "Import", 0));
        menuExport->setTitle(QApplication::translate("parVIEW", "Export", 0));
        secToolBar->setWindowTitle(QApplication::translate("parVIEW", "toolBar", 0));
    } // retranslateUi

};

namespace Ui {
    class parVIEW: public Ui_parVIEW {};
} // namespace Ui

QT_END_NAMESPACE

#endif // UI_PARVIEW_H
