#include "mainwindow.h"
#include "ui_mainwindow.h"
#include <math.h>
#include <iostream>
#include <algorithm>
double norm(std::vector<double> one, std::vector<double> two)
{
    double res = 0;
    for (auto c1 = one.cbegin(); c1 != one.cend(); c1++)
    {
        for (auto c2 = two.cbegin(); c1 != two.cend(); c1++)
        {
            res += pow(((*c1) * (*c1) - (*c2) * (*c2)), 2);
        }
    }
    return sqrt(res);
}

double f_right(double x, double y)
{
    return fabs(pow(sin(M_PI * x * y), 3));
}
double myu1(double y)
{
    return -pow(y, 2) + 1;
}
double myu2(double y)
{
    return -pow(y, 2) + 1;
}
double myu3(double x)
{
    return fabs(sin(M_PI * x));
}
double myu4(double x)
{
    return fabs(sin(M_PI * x));
}

///////
double myu1_test(double y)
{
    return exp(- pow(y, 2));
}
double myu2_test(double y)
{
    return exp(- pow(y, 2));
}
double myu3_test(double x)
{
    return exp(- pow(x, 2));
}
double myu4_test(double x)
{
    return exp(- pow(x, 2));
}

double acc_u(double x, double y)
{
    return exp(1 - pow(x, 2) - pow(y, 2));
}
double f_test(double x, double y)
{
    return 4 * (pow(x, 2) + pow(y, 2) - 1) * exp(1 - pow(x, 2) - pow(y, 2));
}
std::vector<std::vector<double>> solve_Zeidel(std::vector<std::vector<double>> &startSolution,
                                 std::vector<std::vector<double>> &f, unsigned int n, unsigned int m,
                                 unsigned long N, double eps)
{
    unsigned int S = 0;
    double epsMax = 0.0;
    double epsCur = 0.0;
    double a2, k2, h2;
    // std::vector<std::vector<double>> v(n+1, std::vector<double> (m+1));
    double a = -1, b = 1, c = -1, d = 1;
    int i, j;
    double vOld;
    double vNew;
    h2 = -(n/(b-a))*(n/(b-a));
    k2 = -(m/(d-c))*(m/(d-c));
    a2 = -2*(h2+k2);
    auto v = startSolution;
    while (true)
    {
        epsMax = 0;
        for (j = 1; j < m; j++) {
            for (i = 1; i < n; i++) {
                vOld = v[i][j];
                vNew = -(h2*(v[i+1][j] + v[i-1][j]) + k2*(v[i][j-1] + v[i][j-1]));
                vNew += -f[i - 1][j - 1];
                vNew /= a2;
                epsCur = fabs(vOld - vNew);
                if (epsCur > epsMax)
                {
                    epsMax = epsCur;
                }
                v[i][j] = vNew;

            }
            S += 1;

        }
        if ((epsMax < eps) || (S >= N))
        {
            break;
        }
    }
    return v;
}

MainWindow::MainWindow(QWidget *parent) :
    QMainWindow(parent),
    ui(new Ui::MainWindow)
{
    ui->setupUi(this);
}

MainWindow::~MainWindow()
{
    delete ui;
}

void MainWindow::on_pushButton_clicked()
{
    // тестовая задача
    // a = 1, b = 2, c = 2, d = 3

    // Utest = sin(pi*x*y)

    unsigned int n = ui->spinBox->value();
    unsigned int m = ui->spinBox_2->value();
    unsigned int limit = ui->spinBox_3->value();
    double eps = ui->textEdit->toPlainText().toDouble();
    double h = double(2) / double(n);
    double k = double(2) / double(m);
    auto startSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));
    auto f = std::vector<std::vector<double>>(n - 1, std::vector<double> (m - 1, 0));
    for (int i = 0; i < n - 1; ++i)
    {
        for (int j = 0; j < m -1; ++j)
        {
            f[i][j] = f_right(-1 + (i + 1) * h, -1 + (j + 1) * k);
        }
    }
    for (int i = 0; i <= m; ++i)
    {
        startSolution[i][0] = myu1_test(-1 + i * k);
    }
    for (int i = 0; i <= m; ++i)
    {
        startSolution[i][n] = myu2_test(-1 + i * k);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[0][i] = myu3_test(-1 + i * h);
    }
    for (int i = 0; i <= n; ++i)
    {
        startSolution[m][i] = myu1_test(-1 + i * h);
    }

    auto solution = solve_Zeidel(startSolution, f, n, m, limit, eps);
    auto accurateSolution = std::vector<std::vector<double>>(n + 1, std::vector<double> (m + 1, 0));

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            accurateSolution[i][j] = acc_u(-1 + i * h, -1 + j * k);
        }
    }

    ui->tableWidget->setRowCount(m + 1);
    ui->tableWidget->setColumnCount(n + 1);

    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            ui->tableWidget->setItem(i, j, new QTableWidgetItem(QString::number(solution[i][j])));
        }
    }
    auto maxDiffSol = 0.0;
    for (int i = 0; i < n + 1; ++i)
    {
        for (int j = 0; j < m + 1; ++j)
        {
            double currentElem = fabs(solution[i][j] - accurateSolution[i][j]);
            if (currentElem > maxDiffSol){
                maxDiffSol = currentElem;
            }
        }
    }
    std::cout << maxDiffSol << std::endl;
}

void MainWindow::on_pushButton_2_clicked()
{
    //основная задача
     // f = -e^(-xy^2)

}
