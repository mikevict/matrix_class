#include<iostream>
#include<algorithm>
#include<iomanip>
#include<stdlib.h>
#include<map>
#include<cstring>

using namespace std;

const double eps=1e-7;
class Matrix
{
public:
    Matrix(int row,int col);//Ĭ�ϳ�ʼֵȫΪ��
    Matrix& operator+=(Matrix& tmp);
    Matrix& operator-=(Matrix& tmp);
    Matrix& operator*=(Matrix& tmp);
    Matrix& operator*=(int x);
    Matrix& operator/(double x);
    Matrix& operator=(const Matrix& op);
    Matrix& operator=(double *);
    Matrix operator+(Matrix &op);
    Matrix operator-(Matrix &op);
    Matrix operator*(Matrix &op);
    friend void solve(Matrix &op1,Matrix &op2,double* ans);
    friend istream& operator>>(istream &is,Matrix& tmp);
    void operator()(int row,int col,double e);
    friend Matrix& inv(Matrix &op);
    friend int Rank(Matrix& op);
    double det(double** a,int n);
    void Initial_Size(Matrix &op);
    bool  Eig();
    Matrix(Matrix& p);
    void Display_Matrix();
    ~Matrix();
    static Matrix& T(const Matrix &tmp);
    Matrix accompany(Matrix &op);
    int Row();
    int Col();
    double **Data();
    void Row_Ladder();
    Matrix quick_mi(Matrix tmp,int k);
    Matrix operator^(int n);
    void eye();

private:
    int row,col;
    double** data;
};

void Matrix::eye()
{
	for(int i=0;i<row;i++)
	{
		for(int j=0;j<col;j++)
		{
			if(i==j)
			{
				(*this).data[i][j]=1;
			}
		}
	}
}

 

Matrix Matrix::operator^(int n)
{
	Matrix ans(row, col), base(*this);
	ans.eye();
	base.Display_Matrix();
	while(n)
	{
		if(n&1)
			ans *= base;
		n>>=1;
		base *= base;
	}
	return ans;
}


Matrix quick_mi(Matrix tmp,int k)
{
    Matrix x=tmp;
    while(k--)
    {
        x*=tmp;
    }
    return x;
}

void menu()
{
	cout<<"***************************"<<endl; 
    cout<<"1 �ӵȣ�����,�˵�"<<endl;
    cout<<"2 ת��"<<endl;
    cout<<"3 ����"<<endl;
    cout<<"4 ����ʽ"<<endl;
    cout<<"5 �Ӽ���"<<endl;
//  cout<<"6 �ı�ֵ"<<endl;
    cout<<"6 AX=b"<<endl;
    cout<<"7 ����"<<endl;
    cout<<"8 ������ֵ����������"<<endl;
    cout<<"9 �������"<<endl;
    cout<<"***************************"<<endl; 
}
bool  Matrix::Eig()
    {
        if(col!=row)
        {
            cout<<"�����Ƿ���,�޷�������ֵ����������!"<<endl;
            return 0;
        }
        int cnt=10000;//�������� 
        double precision=0.0001;//���þ��� 
        //double *mat=new double[row*col]; 
        Matrix mat(row,col);
        Matrix pvector(row,col);
        
        //double *pvector=new double[row*col];//����������� 
        double *pvalues=new double[row];//�������ֵ 
    //  memset(pvector,0,sizeof(double)*row*col);
    //  cout<<'#'<<endl;
    //  cout<<sizeof(double)*<<endl;
    //  memset(mat,0,sizeof(double)*row*col);
    //  memset(pvalues,0,sizeof(double)*row); 
    //  memcpy(mat,data,sizeof(double)*row*col);
        mat=*this;
    //          for(int i=0;i<row;i++)
    // {
    //  for(int j=0;j<col;j++)
    //  {
    //      cout<<fixed<<setprecision(4)<<data[i*row+j]<<' ';
    //  }
    //  cout<<endl;
    // }
        
//      for(int i=0;i<row;i++)
//      {
//          cout<<pvalues[i]<<' ';
//      }
        for(int i = 0; i < row; i ++) 
        {   
            pvector.data[i][i] = 1.0f; 
            for(int j = 0; j < row; j ++) 
            { 
                if(i != j)   
                    pvector.data[i][j]=0.0f; 
            } 
        }
        

        
        int nCount = 0;     //��������
        while(1)
        {
            //��mat�ķǶԽ������ҵ����Ԫ��
            double dbMax = mat.data[0][1];
            int nRow = 0;
            int nCol = 1;
            for (int i = 0; i < row; i ++)          //��
            {
                for (int j = 0; j < row; j ++)      //��
                {
                    double d = fabs(mat.data[i][j]); 
     
                    if((i!=j) && (d> dbMax)) 
                    { 
                        dbMax = d;   
                        nRow = i;   
                        nCol = j; 
                    } 
                }
            }   
        //  cout<<"@"<<endl; 
            if(dbMax < precision)     //���ȷ���Ҫ�� 
                break;   
            if(nCount > cnt)       //����������������
                break;   
            nCount++;    
            double dbApp = mat.data[nRow][nRow];
            double dbApq = mat.data[nRow][nCol];
            double dbAqq = mat.data[nCol][nCol];     
            //������ת�Ƕ�
            double dbAngle = 0.5*atan2(-2*dbApq,dbAqq-dbApp);
            double dbSinTheta = sin(dbAngle);
            double dbCosTheta = cos(dbAngle);
            double dbSin2Theta = sin(2*dbAngle);
            double dbCos2Theta = cos(2*dbAngle);     
            mat.data[nRow][nRow] = dbApp*dbCosTheta*dbCosTheta + 
                dbAqq*dbSinTheta*dbSinTheta + 2*dbApq*dbCosTheta*dbSinTheta;
            mat.data[nCol][nCol] = dbApp*dbSinTheta*dbSinTheta + 
                dbAqq*dbCosTheta*dbCosTheta - 2*dbApq*dbCosTheta*dbSinTheta;
            mat.data[nRow][nCol] = 0.5*(dbAqq-dbApp)*dbSin2Theta + dbApq*dbCos2Theta;
            mat.data[nCol][nRow] = mat.data[nRow][nCol];     
            for(int i = 0; i < row; i ++) 
            { 
                if((i!=nCol) && (i!=nRow)) 
                { 
                    //int u = i*row + nRow; //p  
                    //int w = i*row + nCol; //q
                    dbMax = mat.data[i][nRow]; 
                    mat.data[i][nRow]= mat.data[i][nCol]*dbSinTheta + dbMax*dbCosTheta; 
                    mat.data[i][nCol]= mat.data[i][nCol]*dbCosTheta - dbMax*dbSinTheta; 
                } 
            }    
            for (int j = 0; j < row; j ++)
            {
                if((j!=nCol) && (j!=nRow)) 
                { 
                    //int u = nRow*row + j; //p
                    //int w = nCol*row + j; //q
                    dbMax = mat.data[nRow][j]; 
                    mat.data[nRow][j]= mat.data[nCol][j]*dbSinTheta + dbMax*dbCosTheta; 
                    mat.data[nCol][j]= mat.data[nCol][j]*dbCosTheta - dbMax*dbSinTheta; 
                } 
            }    
            //������������
            for(int i = 0; i < row; i ++) 
            { 
            //  int u = i*row + nRow;   //[i][nRow] //p   
            //  int w = i*row + nCol;   //[i][nCol] //q
                dbMax = pvector.data[i][nRow]; 
                pvector.data[i][nRow] = pvector.data[i][nCol]*dbSinTheta + dbMax*dbCosTheta; 
                pvector.data[i][nCol] = pvector.data[i][nCol]*dbCosTheta - dbMax*dbSinTheta; 
            }
     
        }    
        //������ֵ���������Լ�����������������,����ֵ��mat���Խ����ϵ�Ԫ��
//      map<double,int> mapEigen;
        cout<<"��������ֵΪ:";
        for(int i=0;i<row;i++)
        {
            cout<<mat.data[i][i]<<' ';
        //  cout<<endl;
        }
        cout<<endl<<endl<<"������������Ϊ(���д��):"<<endl;
        for(int i=0;i<row;i++)
        {
            for(int j=0;j<col;j++)
            {
                printf("%.4lf    ",pvector.data[i][j]);
            }
            cout<<endl;
        }
        return 1;
    }






Matrix& Matrix::operator=(double *tmp)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]=*(tmp+col*i+j);
        }
    }
    return *this;
}
//a is data
void Matrix::Row_Ladder()
{
    for(int i=0;i<col;i++)
    {
        int k;
        for(k=i;k<row;k++)
        {
            if(abs(data[k][i])>eps)
            {
                break;
            }
        }
        if(k<row)
        {
            double temp;
            for(int j=i;j<col;j++)
            {
//              temp=data[i][j];
//              data[i][j]=data[k][j];
//              data[k][j]=temp;
                swap(data[i][j],data[k][j]);
            }
            for(int j=i+1;j<row;j++)
            {
                double t=data[j][i]/data[i][i];
                for(int k=i;k<col;k++)
                {
                    data[j][k]=data[j][k]-t*data[i][k];
                }

            }
        }
    }
}
int Rank(Matrix& op)
{
    Matrix tmp(op);
    tmp.Row_Ladder();
//  tmp.Display_Matrix();
    int count=0;
    for(int i=0;i<tmp.row;i++)
    {
        for(int j=0;j<tmp.col;j++)
        {
            if(abs(tmp.data[i][j])>eps)
            {
                count++;
                break;
            }
        }
    }
    return count;
}



 void solve(Matrix &op1,Matrix &op2,double *ans)
{
    //����ϲ�
    if(op1.row!=op2.row)
    {
        cout<<"solve Ax=b error"<<endl;
        exit(0);
    }
    int y=op1.col+op2.col+3;
    int x=op1.row+3;
    int N=op1.col;
    double A[x][y];
    for(int i=0;i<op1.row;i++)
        for(int j=0;j<op1.col;j++)
        {
            A[i+1][j+1]=op1.data[i][j];
        }
        
        for(int j=op1.col,k=0;j<op1.col+op2.col && k<op2.col;j++,k++)
        {
            for(int i=0;i<op1.row;i++)
            {
                A[i+1][j+1]=op2.data[i][k];
            }
        }
    for(int i=1;i<=N;i++)
    {
        if(abs(A[i][i]-0)<eps)
        {
            cout<<"No Solution"<<endl;
            exit(0);
        }
        double tmp=A[i][i];
        for(int j=i;j<=N+1;j++)
        {
            A[i][j]/=tmp;
        }
        for(int j=i+1;j<=N;j++)
        {
            tmp=A[j][i];
            for(int k=i;k<=N+1;k++)
            {
                A[j][k]-=A[i][k]*tmp;
            }
        }
    }
    for(int i=N;i>=1;i--)
    {
        ans[i]=A[i][N+1];
        for(int j=N;j>=i+1;j--)
        {
            ans[i]-=ans[j]*A[i][j];
        }
    }
}



void Matrix::operator()(int t_row, int t_col,double e)
{
    (*this).Data()[t_row][t_col]=e;
    return;
}

Matrix& inv(Matrix &op)
{

    Matrix temp(op.Row(),op.Col());
    temp=op.accompany(op);
    op=temp/(op.det(op.Data(),op.Row()));
    return op;
}


Matrix Matrix::accompany(Matrix &op)
{
    int n=op.row;
     Matrix tmp(n,n);
    Matrix temp(n,n);
    for(int i=0;i<n;i++)
    {
        for(int j=0;j<n;j++)
        {
            for(int k=0;k<n-1;k++)
            {
                for(int t=0;t<n-1;t++)
                {
                    temp.data[k][t]=op.data[k>=i?k+1:k][t>=j?t+1:t];
                }
            }//������󣬶�Ӧ����
            tmp.data[j][i]=det(temp.data,n-1)*(((i+j+2)%2)?-1:1);
        }
    }
    return tmp;
}




double** Matrix::Data()
{
    return data;
}
int Matrix::Col()
{
    return col;
}

int Matrix::Row()
{
    return row;
}

void Matrix::Display_Matrix()
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            if(data[i][j]-int(data[i][j])==0)
            cout<<data[i][j]<<' ';
            else
            cout<<fixed<<setprecision(4)<<data[i][j]<<' ';
        }
        cout<<endl;
    }
    return;
}


double Matrix::det(double **D,int n)
{
    double d=0;
    if(n==1)
    {
        d=D[0][0];
    }
    else
        if(n==2)
        {
            d=D[0][0]*D[1][1]-D[0][1]*D[1][0];
        }
    else
    {
        for(int i=0;i<n;i++)
        {
            double **M;
            M=new double*[n-1];
            for(int j=0;j<n-1;j++)
                M[j]=new double[n-1];

            for(int j=0;j<n-1;j++)
            {
                for(int k=0;k<n-1;k++)
                {
                    M[j][k]=D[j+1][k>=i?k+1:k];
                }
            }

            if(D[0][i])
            {
                d+=D[0][i]*det(M,n-1)*(((i+2)%2)?-1:1);
            }
            for(int i=0;i<n-1;i++)
            {
                free(M[i]);
            }
                free(M);
        }
        
    }
    //cout<<d<<endl;
    return d;
}


Matrix Matrix::operator+(Matrix &op)
{
    //cout<<"+ ing"<<endl;
    if(col!=op.col||row!=op.row)
    {
        cout<<"+ error"<<endl;
        exit(0);
    }
     Matrix res(row,col);
    Initial_Size(res);
    for(int i=0;i<op.row;i++)
    {
        for(int j=0;j<op.col;j++)
        {
                res.data[i][j]=data[i][j]+op.data[i][j];
        }
    }
    return res;
}

Matrix Matrix::operator-(Matrix &op)
{
    if(col!=op.col||row!=op.row)
    {
        cout<<"- error"<<endl;
        exit(0);
    }
     Matrix res(row,col);
    Initial_Size(res);
    for(int i=0;i<op.row;i++)
    {
        for(int j=0;j<op.col;j++)
        {
                res.data[i][j]=data[i][j]-op.data[i][j];
        }
    }
    return res;
}
Matrix Matrix::operator*(Matrix &op)
{
    if(col!=op.row)
    {
        cout<<"* error"<<endl;
        exit(0);
    }
     Matrix res(row,op.col);
    Initial_Size(res);
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<op.col;j++)
        {
            for(int k=0;k<col;k++)
            {
                res.data[i][j]+=data[i][k]*op.data[k][j];
            }
        }
    }
    return res;
}



Matrix& Matrix::T(const Matrix &op)
{
    int tmp_row=op.col;
    int tmp_col=op.row;
    static Matrix tmp(tmp_row,tmp_col);
    for(int i=0;i<op.row;i++)
        {
            for(int j=0;j<op.col;j++)
            {
                    tmp.data[j][i]=op.data[i][j];
            }
        }
        return tmp;
}

istream& operator>>(istream &is,Matrix& tmp)
{
    cout<<"row "<<tmp.row<<endl<<"col "<<tmp.col<<endl;
    cout<<"please input data"<<endl; 
    for(int i=0;i<tmp.row;i++)
    {
        for(int j=0;j<tmp.col;j++)
        {
            is>>tmp.data[i][j];
        }
    }
    return is;
}


void Matrix::Initial_Size(Matrix &op)
{
    op.data=new double*[op.row];//data ��������double*
    for(int i=0;i<op.row;i++)
    {
        op.data[i]=new double[op.col];
    }
}

Matrix::~Matrix()
{
    for(int i=0;i<row;i++)
    {
        delete []data[i];
    }
    delete []data;
//  cout<<'@'<<endl; 
}

Matrix& Matrix::operator=(const Matrix &tmp)
{
//  cout<<"= ing"<<endl;
    if(this==&tmp)
    {
        return *this;
    }
    if(row!=tmp.row||col!=tmp.col)
    {
        for(int i=0;i<row;i++)
        {
            delete []data[i];
        }
        delete []data;
        row=tmp.row;
        col=tmp.col;
        Initial_Size(*this);
    } 

    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]=tmp.data[i][j];
        }
    }
    return *this;
}
Matrix::Matrix(Matrix &p)
{
    row=p.row;
    col=p.col;
//  cout<<"�������캯��ʹ����"<<endl;
        data=new double*[row];//double��ָ���ָ��
        for(int i=0;i<row;i++)
        {
            data[i]=new double[col];
        }
        for(int i=0;i<row;i++)
        {
            for(int j=0;j<col;j++)
            {
                data[i][j]=p.data[i][j];
            }
        }
}

//void Matrix::Display_Matrix()
//{
//  for(int i=0;i<row;i++)
//  {
//      for(int j=0;j<col;j++)
//      {
//          cout<<data[i][j]<<' ';
//      }
//      cout<<endl;
//  }
//  return;
//}


Matrix::Matrix(int row1,int col1):row(row1),col(col1)
{
//  cout<<"���캯��ʹ����"<<endl;
        Initial_Size(*this);
        for(int i=0;i<row;i++)
        {
            for(int j=0;j<col;j++)
            {
                data[i][j]=0;
            }
        }
}


Matrix& Matrix::operator*=(int x)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]*=x;
        }
    }
    return *this;
}
Matrix& Matrix::operator/(double x)
{
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]/=x;
        }
    }
    return *this;
}

Matrix& Matrix::operator*=(Matrix &tmp)//�г��� 
{
    if(col!=tmp.row)
    {
        cout<<"operator*= error"<<endl;
        exit(0);
    }
    double suml=0.0;
    double temp[row][tmp.col];
    for(int j=0;j<row;j++)
    {
        for(int k=0;k<tmp.col;k++)
        {
            double sumline=0.0;
            for(int i=0;i<col;i++)
            {
                sumline+=data[j][i]*(tmp.data[i][k]);
            }
            temp[j][k]=sumline;
        }
    }
        col=tmp.col;
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]=temp[i][j];
        }
    }
    return *this;
}
Matrix& Matrix::operator+=(Matrix& tmp)
{
    //cout<<"+= ing"<<endl;
    if(row!=tmp.row||col!=tmp.col)
    {
        cout<<"operator+= error"<<endl;
        exit(0);
    }
//  cout<<"#"<<endl;
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]=data[i][j]+tmp.data[i][j];
        }
    }
    return *this;
}

Matrix& Matrix::operator-=(Matrix& tmp)
{
    //cout<<"-=ing"<<endl;
    if(row!=tmp.row||col!=tmp.col)
    {
        cout<<"operator-= error"<<endl;
        exit(0);
    }
    for(int i=0;i<row;i++)
    {
        for(int j=0;j<col;j++)
        {
            data[i][j]=data[i][j]-tmp.data[i][j];
        }
    }
    return *this;
}

int main()
{
//    Matrix text0(3,3);
//    Matrix text1(3,1);
//    Matrix text2(3,3);  
//    cin>>text0;
//    text0.Display_Matrix();
//    text0(1,1,999);
//    puts("");
//    text0.Display_Matrix();
     
    while(1)
{
    menu();
    char op[5];
    cin>>op;
        
        if(op[0]=='1')
        {

            int r,c;
            cout<<"set up two Matrix input row and col num"<<endl;

            cin>>r>>c;
            Matrix t1(r,c),t2(r,c);//�����µľ������

            cout<<"input first Matrix content and second content"<<endl;
            cin>>t1>>t2;//������Ϣ

            string in_op;
            cout<<"input in_op += a,-= b,*= c"<<endl;
            cin>>in_op;
            if(in_op[0]=='a')
            {
                t1+=t2;
            }
            else if(in_op[0]=='b')
            {
                t1-=t2;
            }
            else
            {
                t1*=t2;
            }
            cout<<"ans="<<endl;
            t1.Display_Matrix();
//          t1.~Matrix();
//          t2.~Matrix();
        }
        else if(op[0]=='2')
        {

            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ
            cout<<"ans="<<endl;
            Matrix ans=t1.T(t1);
            ans.Display_Matrix();
        //  t1.~Matrix();
        }
        else if(op[0]=='3')
        {
            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ

            Matrix ans=inv(t1);
            cout<<"ans="<<endl;
            ans.Display_Matrix();
        //  t1.~Matrix();
        }
        else if(op[0]=='4')
        {
            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ
            double ans=t1.det(t1.Data(),t1.Row());
            cout<<"ans="<<endl;
            cout<<ans<<endl;
        //  t1.~Matrix();
        }
        else if(op[0]=='5')
        {

            int r,c;
            cout<<"set up two Matrix input row and col num"<<endl;

            cin>>r>>c;
            Matrix t1(r,c),t2(r,c);//�����µľ������

            cout<<"input first Matrix content and second content"<<endl;
            cin>>t1>>t2;//������Ϣ

            string in_op;
            Matrix ans(r,c);
            cout<<"input in_op + a,- b,* c"<<endl;
            cin>>in_op;
            if(in_op[0]=='a')
            {
                ans=t1+t2;
            }
            else if(in_op[0]=='b')
            {
                ans=t1-t2;
            }
            else
            {
                ans=t1*t2;
            }
            cout<<"ans="<<endl;
            ans.Display_Matrix();
//          t1.~Matrix();
//          t2.~Matrix();
        }
        else if(op[0]=='6')
        {
            int r,c;
            cout<<"set up A Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ

            cout<<"set up b Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t2(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t2;//������Ϣ
            double ans[t1.Col()+3];
            solve(t1,t2,ans);
            cout<<"ans="<<endl;
            for(int i=1;i<=t1.Col();i++)
            {
                cout<<fixed<<setprecision(4)<<ans[i]<<' ';
            }
//          t1.~Matrix();
//          t2.~Matrix();           
        }
        else if(op[0]=='7')
        {
            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ

            int ans=Rank(t1);
            cout<<"ans="<<endl;
            cout<<ans<<endl;
        //  t1.~Matrix();
        }
        else if(op[0]=='8')
        {
            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);//�����µľ������

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ
            bool a=t1.Eig();
        }
        else if(op[0]=='E')
        {
            break;
        }
        else if(op[0]=='9')
        {
            int r,c;
            cout<<"set up a Matrix input row and col num"<<endl;
            cin>>r>>c;
            Matrix t1(r,c);

            cout<<"input the Matrix content "<<endl;
            cin>>t1;//������Ϣ

            int k;
            cout<<"index"<<endl;
            cin>>k;
            Matrix t2(r,c);
            //t2=quick_mi(t1,k);
            t2=t1^k;
            t2.Display_Matrix();
        }
        system("pause");
    system("cls");
}
//  cout<<'#'<<endl;
    
return 0;   
}

//  cout<<'#'<<endl;
//  text0+=text2;
//  text0.Display_Matrix();
//  text0-=text2;
//  text0.Display_Matrix();
//  Matrix text2;
//  cin>>text0;
//  int temp=Rank(text0);
//  cout<<temp<<endl;
//  (text1+text0+text2).Display_Matrix();
//  text0=text1+text0;
//  
//  solve(text0,text1,ans);
//  for(int i=1;i<=text0.Col();i++)
//  {
//      cout<<fixed<<setprecision(2)<<ans[i]<<' ';
//  }
//  Matrix text1(3,2);
//  cin>>text1;
//  Matrix text3(2,2);
//  text3=text0*text1;
// text3=Matrix::T(text0); 
//  cout<<'#'<<endl;
//  double e;
//  cin>>e;
//  text0(1,1,e);
//  cout<<text0.det(text0.Data(),text0.Row());
//  text0.Display_Matrix(); 
//  text0.Display_Matrix();
//  return 0;
//}


