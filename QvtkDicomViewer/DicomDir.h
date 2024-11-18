#pragma once

#include <QWidget>
#include "ui_DicomDir.h"
#include "DicomDataBase.h"

/*
 * 这个类目前是关于DCMTK的各种测试代码的所在地
 * 重构之后将会成为读取DICOMDIR文件并将DicomDateBase对象进行实例化和赋值
 * 在第一次重构的时候,所有的完整性验证都在这个类完成
 * 在第二次重构的时候,用属性封装所有的DicomDatabase属性,并将完整性验证移动到属性里面
 * 在第三次重构的时候,把界面和业务逻辑分离
 */
class PatientMsg
{
public:
	QString FileID;
	QString PatientID;
	QString PatientName;
	QString BirthDate;
	QString Gender;
public:
	PatientMsg()
	{
		FileID = "NULL";
		PatientID = "NULL";
		PatientName = "NULL";
		BirthDate = "NULL";
		Gender = "NULL";
	}
};
class DicomDir : public QWidget
{
	Q_OBJECT

public:
	DicomDir(QWidget *parent = Q_NULLPTR);
	void InitDirExplorerFromDirPath(QString DicomDirFidlePath);//从dirfile绝对路径初始化
	void InitDirExplorerFromSingleFilePath(QString ImageFilePath);//从单图片文件的绝对路径初始化
	void InitDirExplorerFromSeriesPath(QString SeriesPath);//从series的文件夹路径初始化
	~DicomDir();
private:
	Ui::DicomDir ui;
	DicomDataBase* m_database;
	void ConstructsTable();//构造表格
signals:
	void sendData(QString Id);//信号,通过该信号发出病人的ID
public slots:
	void OnPushOk();
	void OnPushCancel();
};
