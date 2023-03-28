using System;
using System.Collections.Generic;

using static GlobalMembers;
using System.IO;

public class Matrix 
{
	public double[,] mat;
    public Matrix(Matrix m)
    {
        this.rows = m.rows;
        this.cols = m.cols;
        allocateMatrix();
        for (int i = 0; i < rows; ++i)
        {
            for (int j = 0; j < cols; ++j)
            {
                mat[i,j] = m.mat[i,j];
            }
        }
    }
   /* public Matrix(List<List<double>> listlist)
	{
        int count = 0;
		//rows = (int)(listlist.GetEnumerator);
        foreach(var enumeration in listlist)
        {
            count++;
        }
        rows = count;
        cols = (int)listlist.Count;

		mat = new double [rows][];

		for (int i = 0; i < rows; i++)
		{
			
			for (int j = 0; j < cols; j++)
			{
				mat[i,j] = (listlist.GetEnumerator() ).begin()[j];
			}
		}
	}*/

	public Matrix(int rows, int cols)
	{
		this.rows = rows;
		this.cols = cols;
        allocateMatrix();
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				mat[i,j] = 0;
			}
		}
	}

	public Matrix(double[][] a, int rows, int cols)
	{
		this.rows = rows;
		this.cols = cols;
        allocateMatrix();
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				mat[i,j] = a[i][j];
			}
		}
	}
    

    public Matrix()
	{
		this.rows = 1;
		this.cols = 1;
		allocateMatrix();
		mat[0,0] = 0;
	}


	
	/*public void Dispose()
	{
		for (int i = 0; i < rows; ++i)
		{
			Arrays.DeleteArray(mat[i]);
		}
		Arrays.DeleteArray(mat);
	}
*/
	
   
	public double functorMethod(int x, int y)
	{
		return mat[x,y];
	} // only for debug

	public static Matrix operator + (Matrix m1, Matrix m2)
	{
        Matrix temp = new Matrix(m1.rows,m1.cols);
        if(m1.rows == m2.rows && m2.cols == m1.cols) { 
        for (int i = 0; i < temp.rows; ++i)
		{
			for (int j = 0; j < temp.cols; ++j)
			{
                temp.mat[i,j] += m2.mat[i,j];
			}
		}
		return temp;
        }
        else
        {
            Console.WriteLine("Not same dimensions.");
            return null;
        }
	}

    public static Matrix operator -(Matrix m1, Matrix m2)
    {
        Matrix temp = new Matrix(m1.rows, m1.cols);
        if (m1.rows == m2.rows && m2.cols == m1.cols)
        {
            for (int i = 0; i < temp.rows; ++i)
            {
                for (int j = 0; j < temp.cols; ++j)
                {
                    temp.mat[i,j] -= m2.mat[i,j];
                }
            }
            return temp;
        }
        else
        {
            Console.WriteLine("Not same dimensions.");
            return null;
        }
    }

    
    public static double[,] operator *(Matrix m1, Matrix m2)
    {
		Matrix temp = new Matrix(m1.rows, m1.cols);
		double[,] temp1 = new double[m1.rows,m1.cols];
		for (int i = 0; i < m1.rows; ++i)
		{
			for (int j = 0; j < m2.cols; ++j)
			{
				for (int k = 0; k < m1.cols; ++k)
				{
					temp1[i,j] += 
						(m1.mat[i,k] * 
						m1.mat[k,j]);
				}
			}
		}
		return temp1;
	}

	public static Matrix operator * (Matrix m, double num)
	{
        for (int i = 0; i < m.rows; ++i)
        {
            for (int j = 0; j < m.cols; ++j)
            {
                m.mat[i,j] *= num;
			}
		}
		return m;
	}
	public static double[] matrix_vector_mult(double[] c, Matrix A, double[] b, int n, int m)
	{
		int i;
		int j;
		int k;
		double sum;

		for (i = 0; i < n; i++)
		{
			sum = 0.0;
			for (k = 0; k < m; k++)
			{
				sum = sum + A.mat[i, k] * b[k];
			}
			c[i] = sum;
			
		}
		return c;
	}
    public static double[] operator * (Matrix A, double[] b)
    {
        int i;
        int j;
        int k;
        double sum;
		double[] c = new double[N];

        for (i = 0; i < n; i++)
        {
            sum = 0.0;
            for (k = 0; k < m; k++)
            {
                sum = sum + A.mat[i, k] * b[k];
            }
            c[i] = sum;

        }
        return c;
    }
    public static double[] operator *( double[] b,Matrix A)
    {
        int i;
        int j;
        int k;
        double sum;
        double[] c = new double[N];

        for (i = 0; i < n; i++)
        {
            sum = 0.0;
            for (k = 0; k < m; k++)
            {
                sum = sum + A.mat[i, k] * b[k];
            }
            c[i] = sum;

        }
        return c;
    }
    public static Matrix operator / (Matrix m ,double num)
	{
		for (int i = 0; i < m.rows; ++i)
		{
			for (int j = 0; j < m.cols; ++j)
			{
                m.mat[i,j] /= num;
			}
		}
		return m;
	}


	public Matrix transpose()
	{
		Matrix ret = new Matrix(cols, rows);
		for (int i = 0; i < rows; ++i)
		{
			for (int j = 0; j < cols; ++j)
			{
				//ret.functorMethod(j, i) = mat[i,j];
                ret.setMatrix(j, i, mat[i,j]);

            }
		}
		return ret;
	}

	public void copy_submatrix(Matrix Source, int istart, int depth, int jstart, int length)
	{
		int i;
		int j;
		int k;
		for (i = istart; i < depth; i++)
		{
			for (j = jstart; j < (jstart + length); j++)
			{
				this.mat[i - istart,j - jstart] = Source.mat[i,j];
			} 
		}
	}

	public void erase_epsilons_matrix(int m, int n)
	{
		int i;
		int j;

		for (i = 0; i < m; i++)
		{
			for (j = 0; j < n; j++)
			{
				if (Math.Abs(this.mat[i,j]) < epsilon)
				{
					this.mat[i,j] = 0.0;
				}
			}
		}
	}

	public int find_exiting_id(double[] x, int enter_id, int n, int m)
	{

		int i;
		int j;
		int temp_min_index = 0;
		int init_flag;
		int q=0;
		int unbounded_flag;
		double temp_min = 0;
		double temp;

		for (i = 0; i < n; i++)
		{
			if (d[i] == enter_id)
			{
				q = i;
			}
		}

		init_flag = 0;
		// CHANGE
		unbounded_flag = 1;
		// END OF CHANGE
		for (i = 0; i < m; i++)
		{
            Console.WriteLine("y[" + i + "][" + q + "]" + "=" + mat[i, q] + ",x[ =" + i + "] =" + x[i] + " " + x[i]);
            Console.WriteLine("init_flag = "+ init_flag);
			if (this.mat[i,q] > 0.0)
			{
				// CHANGE
				unbounded_flag = 0;
				// END OF CHANGE
				temp = x[i] / this.mat[i,q];
                Console.WriteLine("i = " + i + ", temp = " + temp);
                if (init_flag == 0)
				{
					temp_min = temp;
					temp_min_index = i;
					init_flag = 1;
				} // if
				else
				{
					if (temp < temp_min)
					{
						temp_min = temp;
						temp_min_index = i;

					} // if
				}

			} // if
            Console.WriteLine("temp_min_index  = " + temp_min_index + ", temp_min  =" + temp_min);

        } // for
        Console.WriteLine("unbounded flag =" + unbounded_flag);

        if (unbounded_flag == 1)
		{
			Console.WriteLine("Unbounded linear program!\n");
            Environment.Exit(0);

        } // if


		
		return temp_min_index;

	} // find_exiting_id

 
    public void inv_gaussian(Matrix A, int n)
    {

        int i;
        int j;
        int k;
        int p;
        double MaxValue;
        double RelativeValue;
        double temp;

        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                W.mat[i, j] = A.mat[i, j];

        for (i = 0; i < n; i++)
            for (j = n; j < 2 * n; j++)
                W.mat[i, j] = 0.0;

        for (i = 0; i < n; i++)
            W.mat[i, n + i] = 1.0;

        Console.WriteLine("Before loop W: ");
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < 2 * n; j++)
                Console.Write(W.mat[i, j] + " ");
            Console.WriteLine();
        } // for


        for (k = 0; k < n; k++)
        {
            Console.WriteLine("k = " + k);
            p = k;
            MaxValue = Math.Abs(W.mat[k, k]);
            for (i = k + 1; i < n; i++)
                if (Math.Abs(W.mat[i, k]) > MaxValue)
                {
                    p = i;
                    MaxValue = Math.Abs(W.mat[i, k]);
                }// if

            Console.WriteLine("p =" + p + " , k = " + k);

            if (p != k)
            {
                swap_rows(W, n, k, p);
            } // if
            RelativeValue = W.mat[k, k];
            Console.WriteLine("RelativeValue = " + RelativeValue);
            W.mat[k, k] = 1.0;


            for (j = k + 1; j < 2 * n; j++)
            {
                temp = W.mat[k, j] / RelativeValue;
                if (Math.Abs(temp) < epsilon)
                    W.mat[k, j] = 0.0;
                else
                    W.mat[k, j] = temp;
            } // for

            for (i = 0; i < n; i++)
            {
                if (i != k)
                {
                    RelativeValue = W.mat[i, k];
                    W.mat[i, k] = 0.0;
                    for (j = k + 1; j <= 2 * n; j++)
                    {
                        temp = W.mat[i, j] - RelativeValue * W.mat[k, j];
                        if (Math.Abs(temp) < epsilon)
                            W.mat[i, j] = 0.0;
                        else
                            W.mat[i, j] = temp;
                    } // for
                } // if
            } // for


            Console.WriteLine(" W: ");
            for (i = 0; i < n; i++)
            {
                for (j = 0; j < 2 * n; j++)
                    Console.Write(W.mat[i, j] + " ");
                Console.WriteLine();
            } // for


        } /* for */


        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                mat[j, i] = W.mat[j, i + n];
        } // for

        Console.WriteLine("BI:");
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < n; j++)
                Console.Write(mat[i, j] + " ");
            Console.WriteLine();
        } // for


        Console.WriteLine("W:");
        for (i = 0; i < n; i++)
        {
            for (j = 0; j < 2 * n; j++)
                Console.Write(W.mat[i, j] + " ");
            Console.WriteLine();
        } // for




    } /*  gaussian */


    public int getNumOfColumds()
	{
		return cols;
	}

	public int getNumOfRows()
	{
		return rows;
	}


	public void setColumnsNum(int columnsNum)
	{
		cols = columnsNum;
	}

	public void setRowsNum(int rowsNum)
	{
		rows = rowsNum;
	}
    public void setMatrix(int r, int c,double value)
    {
        this.mat[r,c] = value;
    }
	public double[,] getMatPointer()
	{ // only for debug
		return mat;
	}

	private int rows;
	private int cols;

	

	private void allocateMatrix()
	{
		mat = new double[this.rows,this.cols];
		
	}
    public void print()
    {
        int i, j;
            for (i = 0; i < rows; i++)
            for (j = 0; j < cols; j++)
                Console.Write(mat[i, j] + " ");
    }
  
}