#ifndef RESULTTABLE_H
#define RESULTTABLE_H

#include <QWidget>
#include "algebra\vector.hpp"
#include <QModelIndex>
#include "dh_types.h"

QT_BEGIN_NAMESPACE
class QTableWidget;
class QTableWidgetItem;
QT_END_NAMESPACE

class resultTable : public QWidget
{
	Q_OBJECT

public:
	resultTable(QWidget *parent = 0);
	~resultTable();

	void reset();
	void setTable(resultSet& ref_result, algebra::vector<resultSet> &results);
	//void resetTable() { ptable->clearContents(); }
	//void setTable(vector3f* data, std::map<QString, particleObject> *particleObjects);
	private slots:
	void actionClick(QModelIndex);

private:
	QTableWidget *rtable;
};

#endif