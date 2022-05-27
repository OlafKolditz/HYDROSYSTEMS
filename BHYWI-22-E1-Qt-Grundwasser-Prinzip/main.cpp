#include <QApplication>
#include "dialog.h"
#include <QPixmap>
#include <QSplashScreen>
#include <QTime>

int main(int argc, char *argv[])
{
  QApplication a(argc, argv);
  //...........................................................
  QPixmap pixmap ("hydrosystemanalyse.png");
  QSplashScreen splash(pixmap);
  splash.show();
  splash.showMessage(QObject::tr("Übung HSA#1.GW wird geladen..."), Qt::black);
  QTime dieTime = QTime::currentTime().addSecs(3); while( QTime::currentTime() < dieTime )  \
                  QCoreApplication::processEvents(QEventLoop::AllEvents, 1000);
  // Warte-Funktion (http://lists.trolltech.com/qt-interest/2007-01/thread00133-0.html)
  //...........................................................
  Dialog w;
  splash.finish(&w); //test
  w.setWindowTitle("Groundwater Simulator");
  w.setFixedWidth(500);
  w.show();
  //...........................................................
  return a.exec();
}

