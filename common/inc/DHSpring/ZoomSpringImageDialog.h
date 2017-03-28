#ifndef ZOOMSPRINGIMAGEDIALOG_H
#define ZOOMSPRINGIMAGEDIALOG_H

#include <QDialog>

namespace Ui {
	class SpringImageZoomDialog;
}

class ZoomSpringImageDialog : public QDialog
{
	Q_OBJECT

public:
	explicit ZoomSpringImageDialog(QWidget* parent = 0);
	~ZoomSpringImageDialog();

	void setValue();
	//Ui::SpringImageZoomDialog *getUI() { return ui; }

	double f_length;
	double equip_length;
	double init_comp;
	double d;
	double stroke;
	double max_comp;
	double inner_diameter;
	double outer_diameter;
	double middle_diameter;

	Ui::SpringImageZoomDialog *ui;
};

#endif