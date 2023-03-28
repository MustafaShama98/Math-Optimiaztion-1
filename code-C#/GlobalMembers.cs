using System;
using System.Collections.Generic;
using System.Diagnostics;
using System.IO;
using System.Text;
public static class GlobalMembers
{
    public static int M = 64;
    public static int N = 128;
    public static double[] C = new double[N];
    public static double[] c = new double[N];
    public static double[] b = new double[N];
    public static double[] b_aux = new double[N];
    public static double[] cb = new double[N];
    public static double[] cbBI = new double[N];
    public static double[] cbBID = new double[M];
    public static double[] cd = new double[N];
    public static double[] rd = new double[N];
    public static double[] BIb = new double[N];
    public static double epsilon = 0;
    public static int[] d = new int[N];
    public static int[] d_aux = new int[N];
    public static int[] basis = new int[N];
    public static int n;
    public static int m;

    public static int Initial_n;
    public static double[] Initial_cb = new double[N];
    public static double[] Initial_cd = new double[N];
    public static double[] Initial_c = new double[N];
    public static double[] Initial_c_aux = new double[N];
    public static int[] Initial_basis = new int[N];
    public static int[] Initial_d = new int[N];
    public static double[] Initial_BIb = new double[N];
    public static double[] Initial_b = new double[N];
    public static double[] Initial_b_aux = new double[N];
    public static double[] Initial_rd = new double[N];
    public static double[] Initial_cbBI = new double[N];
    public static double[] Initial_cbBID = new double[N];
  
    public static Matrix BID = new Matrix(N, N);
    public static Matrix W = new Matrix(N, N);
    public static Matrix BI = new Matrix(N, N);
    public static Matrix BIA_aux = new Matrix(N, N);
    public static Matrix Initial_D = new Matrix(N, N);
    public static Matrix A_aux = new Matrix(N, N);
    public static Matrix D = new Matrix(N, N);
    public static Matrix A = new Matrix(N, N);
    public static Matrix B = new Matrix(N, N);
    public static Matrix Initial_BID = new Matrix(N, N);
    public static Matrix Initial_BI = new Matrix(N, N);
    public static Matrix Initial_B = new Matrix(N, N);
    public static Matrix Initial_C = new Matrix(N, N);
    public static Matrix Initial_A = new Matrix(N, N);
    public static Matrix Initial_A_aux = new Matrix(N, N);
    public static Matrix Initial_W = new Matrix(N, N);
    public static Matrix Initial_BIA_aux = new Matrix(N, N);

	//A function to delete extra white spaces in text files so it matches the format to parse the input data
    public static string TrimAllExtraWhiteSpaces( string input)
    {
        if (string.IsNullOrEmpty(input))
        {
            return input;
        }

        var current = 0;
        char[] output = new char[input.Length];
        var charArray = input.ToCharArray();

        for (var i = 0; i < charArray.Length; i++)
        {
            if (!char.IsWhiteSpace(charArray[i]))
            {
                if (current > 0 && i > 0 && char.IsWhiteSpace(charArray[i - 1]))
                {
                    output[current++] = ' ';
                }
                output[current++] = charArray[i];
            }
        }

        return new string(output, 0, current);
    }
    static void Main(string[] args)
        {

            int i;
            int j;;
            int n_p_m;
            int itemp;
            string path = "..\\..\\";
		string temp;
        int lineNumber = 0;
		Console.Write("Please insert name of the text file (without .txt): ");
		string fileName = Console.ReadLine();
        fileName=path + fileName +".txt";
        foreach (string line in System.IO.File.ReadLines(fileName)) 
        {
            int index;
            if (lineNumber == 1)
            {
                temp = TrimAllExtraWhiteSpaces(line);
                string[] textPart = temp.Split(' ');
                m = Int32.Parse(textPart[0]);
                n = Int32.Parse(textPart[1]);
                Console.WriteLine("m, n = " + m + ", " + n);
            }

            if (lineNumber == 3)//C
            {
				temp=TrimAllExtraWhiteSpaces(line);
                string[] cStr = temp.Split(' ');
               
                Console.WriteLine();
                Console.Write("C:" + " ");

                for (index = 0; index < n; index++)
                {
                    c[index] = double.Parse(cStr[index]);
                    Console.Write(c[index] + " ");

                }
            }
            if (lineNumber >= 5 && lineNumber <= 4 + m)//A
            {
                temp = TrimAllExtraWhiteSpaces(line);
                string[] AStr = temp.Split(' ');
                Console.WriteLine();
                Console.Write("A:" + " ");
                for (index = 0; index < n; index++)
                {
                   // AStr.spli
                    A.mat[lineNumber - 5, index] = double.Parse(AStr[index]);
                    Console.Write(double.Parse(AStr[index]) + " ");

                }
            }
            if (lineNumber == 6 + m)//b
            {
                temp = TrimAllExtraWhiteSpaces(line);
                string[] bStr = temp.Split(' ');
                Console.WriteLine();
                Console.Write("b:" + " ");

                for (index = 0; index < m; index++)//size m
                {
                    b[index] = double.Parse(bStr[index]);
                    Console.Write(b[index] + " ");

                }
            }
            if (lineNumber == 8 + m)//epsilon
            {
                temp = TrimAllExtraWhiteSpaces(line);
                string[] epStr = temp.Split(' ');
                Console.WriteLine();
                epsilon = double.Parse(epStr[0]);

            }
            lineNumber++;
        }

        Initial_n = n + m;
        n_p_m = n + m;

        copy_matrix(A_aux, A, n, m);

        Console.WriteLine("epsilon = " + epsilon);

        Console.WriteLine(" A: ");
        print_original_system();


        copy_to_initial_matrix();

        Console.WriteLine("Initial_A:");
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n_p_m; j++)
                Console.Write(Initial_A.mat[i, j]+" " );
            Console.WriteLine();
        } // for



        for (i = 0; i < m; i++)
            Initial_basis[i] = i + n;
        for (i = 0; i < n; i++)
            Initial_c[i] = 0.0;
        for (i = n; i < Initial_n; i++)
            Initial_c[i] = 1.0;
        for (i = 0; i < m; i++)
            Initial_b[i] = b[i];
        for (i = 0; i < m; i++)
            Initial_b_aux[i] = b[i];

        Console.WriteLine("Initial_basis:");
        for (i = 0; i < m; i++)
            Console.Write(Initial_basis[i]);
        Console.WriteLine();

        Console.WriteLine("Initial_c:");
        for (i = 0; i < Initial_n; i++)
            Console.Write(Initial_c[i]);
        Console.WriteLine();


        Console.WriteLine("Initial_b:");
        for (i = 0; i < m; i++)
            Console.Write(Initial_b[i]);
        Console.WriteLine();


        Initial_simplex_algorithm();


        for (i = 0; i < m; i++)
        {
            itemp = Initial_basis[i];
            basis[i] = itemp;
            if (itemp >= n)
            {
                print_no_solution();
                Console.WriteLine("Press any key to continue . . .");
                Console.ReadKey();
            } 

        }
        print_initial_solution();
            simplex_algorithm();
            bublesort_d(basis, BIb, m);

            print_solution();
        Console.WriteLine("Press any key to continue . . .");
        Console.ReadKey();
    } // main

    

    public static void bublesort(int[] arr, int n)
	{
	  int i;
	  int j;
	  int limit;
	  int flag;
	  int temp;

	  flag = 1;
	  for (i = 0;(i < n) && (flag == 1); i++)
	  {
		flag = 0;
		limit = n - i - 1;
			for (j = 0; j < limit; j++)
			{
			   if (arr[j] > arr[j + 1])
			   {
					 flag = 1;
					 temp = arr[j];
					 arr[j] = arr[j + 1];
					 arr[j + 1] = temp;
			   } // if
			}
	  } // for

	} // bublesort
	public static void bublesort_d(int[] arr, double[] darr, int n)
	{
	  int i;
	  int j;
	  int limit;
	  int flag;
	  int temp;
	  double dtemp;

	  flag = 1;
	  for (i = 0;(i < n) && (flag == 1); i++)
	  {
		flag = 0;
		limit = n - i - 1;
			for (j = 0; j < limit; j++)
			{
			   if (arr[j] > arr[j + 1])
			   {
					 flag = 1;
					 dtemp = darr[j];
					 darr[j] = darr[j + 1];
					 darr[j + 1] = dtemp;
					 temp = arr[j];
					 arr[j] = arr[j + 1];
					 arr[j + 1] = temp;
			   } // if
			}
	  } // for

	} // bublesort

	public static void compute_cb_cd()
	{

	  int i;



	  for (i = 0; i < m; i++)
	  {
	   cb[i] = c[d[i]];
	  }

	  for (i = m; i < n; i++)
	  {
	   cd[i - m] = c[d[i]];
	  }

	  Console.Write("d:\n");
	  for (i = 0; i < n; i++)
	  {
	   Console.Write(d[i]+" ");
	  }
	  Console.Write("\n");

	  Console.Write("cb:\n");
	  for (i = 0; i < m; i++)
	  {
	   Console.Write(cb[i]+ " ");
	  }
	  Console.Write("\n");

	  Console.Write("cd:\n");
	  for (i = 0; i < n; i++)
	  {
	   Console.Write(cd[i]+ " ");
	  }
	  Console.Write("\n");


	} // compute_cb_cd()



	public static void compute_Initial_cb_Initial_cd()
	{

	  int i;



	  for (i = 0; i < m; i++)
	  {
		  Initial_cb[i] = Initial_c[Initial_d[i]];
	  }

	  for (i = m; i < Initial_n; i++)
	  {
	   Initial_cd[i - m] = Initial_c[Initial_d[i]];
	  }

	  Console.Write("Initial_d:\n");
	  for (i = 0; i < Initial_n; i++)
	  {
	   Console.Write( Initial_d[i]+ " ");
	  }
	  Console.Write("\n");

	  Console.Write("Initial_cb:\n");
	  for (i = 0; i < m; i++)
	  {
	   Console.Write(Initial_cb[i]+ " ");
	  }
	  Console.Write("\n");

	  Console.Write("Initial_cd:\n");
	  for (i = 0; i < (Initial_n - m); i++)
	  {
	   Console.Write( Initial_cd[i]+ " ");
	  }
	  Console.Write("\n");


	} // compute_cb_cd()



	public static void copy_matrix(Matrix Dest, Matrix Source, int n, int m)
	{
	 int i;
	 int j;

	  for (i = 0; i < m; i++)
	  {
	   for (j = 0; j < n; j++)
	   {
		  Dest.mat[i,j] = Source.mat[i,j];
	   }
	  }


	} // copy_matrix


	public static void copy_submatrix(double[][] Dest, double[][] Source, int istart, int depth, int jstart, int length)
	{
	  int i;
	  int j;
	  int k;

	  for (i = istart; i < depth; i++)
	  {
		for (j = jstart; j < (jstart + length); j++)
		{
		  Dest[i - istart][j - jstart] = Source[i][j];
		}
	  }


	} // copy_submatrix



	public static void copy_to_initial_matrix()
	{
	  int i;
	  int j;
	  int k;

	  for (i = 0; i < m; i++)
	  {
		for (j = 0; j < n; j++)
		{
				Initial_A.mat[i,j] = Initial_A_aux.mat[i, j] = A.mat[i, j];
		}
	  }
	//      Initial_A[i][j] = Initial_A_aux[i][j] = A[i][j-m];

	  for (i = 0; i < m; i++)
	  {
		for (j = n; j < n + m; j++)
		{
		 if (i == (j - n))
		 {
			Initial_A.mat[i,j] = Initial_A_aux.mat[i,j] = 1.0;
		 }
		 else
		 {
			Initial_A.mat[i,j] = Initial_A_aux.mat[i,j] = 0.0;
		 }
		}
	  }
	  Initial_A.setColumnsNum(n + m);
	  Initial_A.setRowsNum(m);
	  Initial_A_aux.setColumnsNum(n + m);
	  Initial_A_aux.setRowsNum(m);
	} // copy_to_initial_matrix

	public static void copy_vector(double[] Dest, double[] Source, int n)
	{
	 int i;
	 int j;

	  for (i = 0; i < n; i++)
	  {
		  Dest[i] = Source[i];
	  }


	} // copy_vector



	public static void erase_epsilons_matrix(double[,] dmat, int m, int n)
	{
	  int i;
	  int j;

	  for (i = 0; i < m; i++)
	  {
		for (j = 0; j < n; j++)
		{
		   if (Math.Abs(dmat[i,j]) < epsilon)
		   {
			   dmat[i,j] = 0.0;
		   }
		}
	  }


	} // erase_epsilons_matrix

	



	public static void erase_epsilons_vector(double[] darray, int n)
	{
	  int i;

	  for (i = 0; i < n; i++)
	  {
		if (Math.Abs(darray[i]) < epsilon)
		{
			darray[i] = 0.0;
		}
	  }


	} // erase_epsilons_vector

	



	public static void find_all_negative_rds(int[] neg_ids, ref int no_of_rds)
	{

	  int i;
	  int index;
	  index = 0;
	  for (i = 0; i < (n - m); i++)
	  {
		if (rd[i] < 0)
		{
		   neg_ids[index] = d[m + i];
		   index++;
		   Console.Write("\nXX:i = "+i+" rd["+i+"] = "+rd+", d["+i+"] =  "+d+", index = "+index+"\n");
		} // if
	  }

	  no_of_rds = index;
	} // find_all_negative_rds

	
	public static int find_Initial_exiting_id(Matrix y, double[] x, int enter_id, int n, int m)
	{

	  int i;
	  int j;
	  int temp_min_index = 0;
	  int init_flag;
	  int q=0;
	  double temp_min = 0;
	  double temp;

	  for (i = 0; i < Initial_n; i++)
	  {
		if (Initial_d[i] == enter_id)
		{
			q = i;
		}
	  }

	  init_flag = 0;
	  for (i = 0; i < m; i++)
	  {
	   Console.Write("y["+i+"]["+q+"] = "+ y.mat[i, q] + ", x["+ i + "] = "+ x[i] + "\n");
	   Console.Write("init_flag = "+init_flag+"\n");
	   if (y.functorMethod(i,q) > 0.0)
	   {
		temp = x[i] / y.mat[i, q];
		Console.Write("i = "+i+", temp ="+temp+"\n");
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
		Console.Write("temp_min_index  = "+temp_min_index+", temp_min  = "+temp_min+"\n");

	  } // for

	return temp_min_index;

	} // find_exiting_id

	public static int find_Initial_most_negative()
	{
	// Assumptions: d[i] is original index, rd orderred by i=0,.., n-1 
	//                                                  d[0] ... d[n-1]

	  int most_index;
	  int i;
	  double temp_value;

	  most_index = Initial_d[m];
	  temp_value = Initial_rd[0];
	  for (i = 0; i < (Initial_n - m); i++)
	  {
		if (Initial_rd[i] < temp_value)
		{
		   most_index = Initial_d[m + i];
		   temp_value = Initial_rd[i];
		} // if
	  }

	 return most_index;
	} // find_most_negative


	public static double find_min_value(double[] rd, int n)
	{
	// Assumptions: d[i] is original index, rd orderred by i=0,.., n-1 
	//                                                  d[0] ... d[n-1]

	  int most_index;
	  int i;
	  double temp_value;

	  most_index = d[0];
	  temp_value = rd[0];
	  for (i = 0; i < n; i++)
	  {
		if (rd[i] < temp_value)
		{
		   most_index = d[i];
		   temp_value = rd[i];
		} // if
	  }

	 return temp_value;
	} // find_min_value



	public static int find_most_negative()
	{
	// Assumptions: d[i] is original index, rd orderred by i=0,.., n-1 
	//                                                  d[0] ... d[n-1]

	  int most_index;
	  int i;
	  double temp_value;

	  most_index = d[m];
	  temp_value = rd[0];
	  for (i = 0; i < (n - m); i++)
	  {
		if (rd[i] < temp_value)
		{
		   most_index = d[m + i];
		   temp_value = rd[i];
		} // if
	  }

	 return most_index;
	} // find_most_negative

	public static void Initial_set_d()
	{

	 int i;
	 int j;
	 int pos;
	 int flag;

	 for (i = 0; i < m; i++)
	 {
	   Initial_d[i] = Initial_basis[i];
	 }

	 pos = m;
	 for (i = 0; i < Initial_n; i++)
	 {
	  flag = 1;
	  for (j = 0; (j < m) && (flag == 1); j++)
	  {
	   if (i == Initial_basis[j])
	   {
		 flag = 0;
	   }
	  }
	  if (flag == 1)
	  {
		Initial_d[pos++] = i;
	  }
	 } // for



	} // Initial_set_d
	public static void Initial_simplex_algorithm()
	{
		int i;
		int j;
		int k;
		int optimal_flag;
		int enter_id;
		int exiting_id;
		int itemp;
		int basis_i;
		double dtemp;
		double min_value;
		int count = 0;

		optimal_flag = 0;

		Console.Write(" m = " +m + "Initial_n = " +  Initial_n);

		Console.Write("\nInitial_basis:\n");
		for (i = 0; i < m; i++)
		{
			Console.Write(Initial_basis[i]+ " ");
		}
		Console.Write("\n");


		Console.Write("Initial_A:\n");
		for (i = 0; i < m; i++)
		{
			for (j = 0; j < Initial_n; j++)
			Console.Write(Initial_A.mat[i, j] + " ");
			
		} // for
	
		Console.Write("\n");



		while (optimal_flag == 0)
		{


			bublesort(Initial_basis, m);
			

			Console.Write("\nInitial_basis:");
			for (i = 0; i < m; i++)
			{
				Console.Write(Initial_basis[i] + " ");
			}
			Console.Write("\n");


			Initial_set_d();


			Console.Write("Initial_d:");
			for (i = 0; i < Initial_n; i++)
			{
				Console.Write(Initial_d[i] + " ");
			}
			Console.Write("\n");


			set_Initial_A_aux();

			Console.Write("Initial_A_aux (B, D):\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < Initial_n; j++)
				{
					Console.Write( Initial_A_aux.mat[i,j]+ " ");
				}
				Console.Write(" \n");
			} // for


		
			Initial_B.copy_submatrix(Initial_A_aux, 0, m, 0, m);
			//Initial_B.setColumnsNum(m); // set the used size of the matrix
			//Initial_B.setRowsNum(m); // set the used size of the matrix

			Console.Write("\nInitial_B:\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < m; j++)
					Console.Write( Initial_B.mat[i, j]+ " ");
				Console.Write("\n");
			} 
			
			Console.Write("\n");

			//inv_gaussian(Initial_BI, Initial_B, m); // BI = B-1  
			

            Initial_BI.inv_gaussian(Initial_B, m);
            
            //erase_epsilons_matrix(Initial_BI, m, m);
            Initial_BI.erase_epsilons_matrix(m, m);

			Console.Write("\nInitial_BI:\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < m; j++)
					Console.Write(Initial_BI.mat[i, j] + " ");
				Console.WriteLine();

            } // for  
			Console.Write("\n");
		

			matrix_mult(Initial_BIA_aux, Initial_BI, Initial_A_aux, m, m, Initial_n);
		//	Initial_BIA_aux.mat = Initial_BI * Initial_A_aux;

			Initial_BIA_aux.setColumnsNum(Initial_A_aux.getNumOfColumds()); // set the used size of the matrix
			Initial_BIA_aux.setRowsNum(Initial_BI.getNumOfRows()); // set the used size of the matrix

			//erase_epsilons_matrix(Initial_BIA_aux, m, Initial_n);
			Initial_BIA_aux.erase_epsilons_matrix(m, Initial_n);

			Console.Write("\nInitial_BIA_aux (I, B-1*D):\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < Initial_n; j++)
					Console.Write(Initial_BIA_aux.mat[i, j]+ " ");
				Console.Write("\n");
			} // for  
			

			Console.Write("\nInitial_A_aux (B,D):\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < Initial_n; j++)
					Console.Write( Initial_A_aux.mat[i,j]+ " ");
				Console.Write("\n");
			} // for  
		

			Console.Write("Initial_b:\n");
			for (i = 0; i < m; i++)
			{
				Console.Write(Initial_b[i]+ " ");
			}


			Initial_BIb = Initial_BI * Initial_b;
			erase_epsilons_vector(Initial_BIb, m);

			Console.Write("\nInitial_BIb:\n");

			for (i = 0; i < m; i++)
			{
				Console.Write( Initial_BIb[i]+ " ");
			}
			Console.Write("\n");


		
			Initial_D.copy_submatrix(Initial_A_aux, 0, m, m, Initial_n - m);
			//Initial_D.setColumnsNum(Initial_n - m); // set the used size of the matrix
			//Initial_D.setRowsNum(m); // set the used size of the matrix

			Console.Write("Initial_D:\n");
			for (i = 0; i < m; i++)
			{
				for (j = 0; j < Initial_n - m; j++)
					Console.Write(Initial_D.mat[i,j]+ " ");
				Console.Write("\n");
			} // for  

			compute_Initial_cb_Initial_cd();

			Console.Write("\nInitial_cb:\n");
			for (i = 0; i < m; i++)
			{
				Console.Write( Initial_cb[i]+ " ");
			}
			Console.Write("\n");

			Console.Write("\nInitial_cd:\n");
			for (i = 0; i < (Initial_n - m); i++)
			{
				Console.Write( Initial_cd[i]+ " ");
			}
			Console.Write("\n");


			// cbBI = cb * B-1

			vector_matrix_mult(Initial_cbBI, Initial_cb, Initial_BI, m, m);
			//Initial_cbBI = Initial_cb * Initial_BI;
			erase_epsilons_vector(Initial_cbBI, m);

			Console.Write("\nInitial_cbBI:\n");
			for (i = 0; i < m; i++)
			{
				Console.Write( Initial_cbBI[i]+ " ");
			}
			Console.Write("\n");


			vector_matrix_mult(Initial_cbBID, Initial_cbBI, Initial_D, m, Initial_n - m);
			//Initial_cbBID = Initial_cbBI * Initial_D;
			erase_epsilons_vector(Initial_cbBID, Initial_n - m);


			Console.Write("\nInitial_cbBID:\n");
			for (i = 0; i < (Initial_n - m); i++)
			{
				Console.Write( Initial_cbBID[i] + " ");
			}
			Console.Write("\n");



			vector_subtract(Initial_rd, Initial_cd, Initial_cbBID, Initial_n - m);
			erase_epsilons_vector(Initial_rd, Initial_n - m);

			Console.Write("\nInitial_rd( cd - cbBID ):\n");
			for (i = 0; i < (Initial_n - m); i++)
			{
				Console.Write(Initial_rd[i]+ " ");
			}
			Console.Write("\n\n");


			min_value = find_min_value(Initial_rd, n);
			if (min_value >= 0.0)
			{
				optimal_flag = 1;
			}
			else
			{
				enter_id = find_Initial_most_negative();
				exiting_id = find_Initial_exiting_id(Initial_BIA_aux, Initial_BIb, enter_id, Initial_n, m);
                Console.WriteLine("enter_id  = " + enter_id + ",  exiting_id = " + exiting_id +", Initial_d[exiting_id] = " + Initial_d[exiting_id]);

                Initial_basis[exiting_id] = enter_id;
				Console.Write("\nInitial_basis:\n");
				for (i = 0; i < m; i++)
				{
					Console.Write(Initial_basis[i]+ " ");
				}
				Console.Write("\n");

			} // else

		} // while

	} //Initial_simplex_algorithm


	


	public static void Initial_swap_colums(int i, int j)
	{
	  int k;
	  double temp;



	  for (k = 0; k < m; k++)
	  {
		 //Initial_A_aux[k][i] = Initial_A[k][j];
		 Initial_A_aux.mat[k,i] = Initial_A.mat[k,j];
		 // Initial_A_aux[k][j] = Initial_A[k][i];
		 Initial_A_aux.mat[k,j] = Initial_A.mat[k,i];
	  } // for


	} // swap_colums

	



	public static void inv_gaussian(double[][] B, double[][] A, int n)
	{

	  int i;
	  int j;
	  int k;
	  int p;
	  int itemp;
	  double MaxValue;
	  double RelativeValue;
	  double temp;


	  for (i = 0; i < n; i++)
	  {
		for (j = 0; j < n; j++)
		{
		  W.mat[i,j] = A[i][j];
		}
	  }

	   for (i = 0; i < n; i++)
	   {
		for (j = n; j < 2 * n; j++)
		{
		  W.mat[i,j] = 0.0;
		}
	   }

	   for (i = 0; i < n; i++)
	   {
		 W.mat[i,n + i] = 1.0;
	   }

	   Console.Write("\nBefore loop W: ");
	   for (i = 0; i < n; i++)
	   {
		 for (j = 0; j < 2 * n; j++)
		 {
			 Console.Write(W.mat[i,j]+" ");
		 }
		Console.Write("\n");
	   } // for


	  for (k = 0; k < n; k++)
	  {
	   Console.Write("k ="+  k);
		p = k;
		MaxValue = Math.Abs(W.mat[k,k]);
		for (i = k + 1; i < n; i++)
		{
		 if (Math.Abs(W.mat[i,k]) > MaxValue)
		 {
			 p = i;
			 MaxValue = Math.Abs(W.mat[i,k]);
		 } // if
		}

		 Console.WriteLine("p = "+p, "k = "+k, p, k);

		 if (p != k)
		 {
		   swap_rows(W, n, k, p);
		 } // if
		 RelativeValue = W.mat[k,k];
		Console.WriteLine("RelativeValue = "+ RelativeValue);
		 W.mat[k,k] = 1.0;


		 for (j = k + 1; j < 2 * n; j++)
		 {
			temp = W.mat[k,j] / RelativeValue;
			if (Math.Abs(temp) < epsilon)
			{
			   W.mat[k,j] = 0.0;
			}
			else
			{
				W.mat[k,j] = temp;
			}
		 } // for

		 for (i = 0; i < n; i++)
		 {
		  if (i != k)
		  {
			 RelativeValue = W.mat[i,k];
			 W.mat[i,k] = 0.0;
			 for (j = k + 1; j <= 2 * n; j++)
			 {
			   temp = W.mat[i,j] - RelativeValue * W.mat[k,j];
			   if (Math.Abs(temp) < epsilon)
			   {
				  W.mat[i,j] = 0.0;
			   }
			   else
			   {
				  W.mat[i,j] = temp;
			   }
			 } // for
		  } // if
		 } // for


	   Console.Write(" W: \n");
	   for (i = 0; i < n; i++)
	   {
		 for (j = 0; j < 2 * n; j++)
		 {
			 Console.Write(W.mat[i,j]+" ");
		 }
		Console.Write("\n");
	   } // for


	  } // for


	   for (i = 0; i < n; i++)
	   {
		 for (j = 0; j < n; j++)
		 {
			B[j][i] = W.mat[j,i + n];
		 }
	   } // for

	   Console.Write("\nBI:\n");
	   for (i = 0; i < n; i++)
	   {
		 for (j = 0; j < n; j++)
		 {
			 Console.Write(B[i][j]+ " ");
		 }
		 Console.Write("\n");
	   } // for


	   Console.Write("\nW:\n");
	   for (i = 0; i < n; i++)
	   {
		 for (j = 0; j < 2 * n; j++)
		 {
			 Console.Write(W.mat[i,j]+ " ");
		 }
		 Console.Write("\n");
	   } // for




	} //  gaussian



	public static void matrix_mult(Matrix C, Matrix A, Matrix B, int n, int m, int p)
	{
	  int i;
	  int j;
	  int k;
	  double sum;

	 for (i = 0; i < n; i++)
	 {
	   for (j = 0; j < p; j++)
	   {
		   sum = 0.0;
		   for (k = 0; k < m; k++)
		   {
			  sum = sum + A.mat[i,k] * B.mat[k,j];
		   }
		 C.mat[i,j] = sum;
	   } // for
	 }

	} // matrix_mult



	public static void print_initial_solution()
	{

	  int i;

	  Console.Write("\nInitial basis:\n");
	  for (i = 0; i < m; i++)
	  {
		Console.Write(Initial_basis[i]+ " ");
	  }
	  Console.Write("\n");

	  Console.Write("\nBasic Solution:\n");
	  for (i = 0; i < m; i++)
	  {
            Console.Write(" X" + Initial_basis[i] + "= " + Initial_BIb[i]);
        }
	  Console.Write("\n");

	} // print_initial_solution




	public static void print_no_solution()
	{

	  Console.Write("System A has NO solution\n");

	} // print_no_solution






	public static void print_original_system()
	{
	 int i;
	 int j;

	 Console.Write("Original System:\n");

	 
	 for(i=0; i < m; i++)
	 {
	   for(j=0; j < n; j++)
	     Console.Write(A.mat[i, j] + " ");
	   Console.Write("\n");
	 } 

	} // print_original_system


	public static void print_result()
	{
	 int i;
	 int j;
	 int k;
	 double sum;

	 Console.Write("Optimal Basis:\n");
	 for (i = 0; i < m; i++)
	 {
	   Console.Write( basis[i]+ " ");
	 }
	 Console.Write("\n");


	 Console.Write("Optimal Solution:\n");
	 for (i = 0; i < m; i++)
	 {
            Console.Write(" X" + basis[i] + " = " + BIb[i]);
        }
	 Console.Write("\n");



	} // print_result



//
	internal static int print_simplex_params_count = 0;

	public static void print_simplex_params(double[][] A, double[][] A_aux, double[] c, double[] b, int n, int m, double[][] B, double[][] BID, double[][] D, int[] basis, int[] d, double[] cb, double[] cd)
	{
	  int i;
	  int j;
	//  static int count = 0;

	  Console.WriteLine(" m = "+m+" n = "+n);

	  Console.Write("A:\n");
	  for (i = 0; i < m; i++)
	  {
		 for (j = 0; j < n; j++)
		 {
		   Console.Write(A[i][j]+ " ");
		 }
		 Console.Write("\n");
	  } // for


	  Console.Write("c:\n");
	  for (i = 0; i < n; i++)
	  {
		 Console.Write( c[i]+ " ");
	  }
	  Console.Write("\n");


	  Console.Write("b:\n");
	  for (i = 0; i < m; i++)
	  {
		 Console.Write( b[i]+ " ");
	  }
	  Console.Write("\n");



	  Console.Write("A_aux:\n");
	  for (i = 0; i < m; i++)
	  {
		 for (j = 0; j < n; j++)
		 {
		   Console.Write( A_aux[i][j]+ " ");
		 }
		 Console.Write("\n");
	  } // for




	  Console.Write("B:\n");
	  for (i = 0; i < m; i++)
	  {
		 for (j = 0; j < m; j++)
		 {
		   Console.Write( B[i][j]+ " ");
		 }
		 Console.Write("\n");
	  } // for

	  Console.Write("basis:\n");
	  for (i = 0; i < m; i++)
	  {
		 Console.Write( basis[i]+" ");
	  }
	  Console.Write("\n");

	  print_simplex_params_count++;


	  if (print_simplex_params_count >= 8)
	  {
            Environment.Exit(0);
	  }

	} // print_simplex_params




	public static void print_solution()
	{

	  int i;
	  double sum;
	  double temp;

	  Console.Write("\nbasis:\n");
	  for (i = 0; i < m; i++)
	  {
		Console.Write(basis[i]+ " ");
	  }
	  Console.Write("\n");

	  Console.Write("\nBasic Solution:\n");
	  for (i = 0; i < m; i++)
	  {
            Console.Write(" X" + (basis[i] + 1) + " = " + BIb[i]);
        }
	  Console.Write("\n");

	  Console.Write("\nSolution value:\n");

	  temp = c[basis[0]] * BIb[0];
	  sum = temp;
	  Console.Write(c[basis[0]]+ " *  "+ BIb[0]);

	  for (i = 1; i < m; i++)
	  {
		 temp = c[basis[i]] * BIb[i];
		 sum = sum + temp;
            Console.Write(" + " + c[basis[i]] + " * " + BIb[i]);
             
        } // for

	  Console.WriteLine(" = "+ sum);

	} // print_solution

	






	public static void set_A_aux()
	{

	  int i;
	  int j;
	  int k;


	  for (i = 0; i < n; i++)
	  {
		k = d[i];
		for (j = 0; j < m; j++)
		{
		 A_aux.mat[j,i] = A.mat[j,k];
		} // for 2
	  } // for1



	} // swap_colums


	public static void set_d()
	{

	 int i;
	 int j;
	 int pos;
	 int flag;

	 for (i = 0; i < m; i++)
	 {
	   d[i] = basis[i];
	 }

	 pos = m;
	 for (i = 0; i < n; i++)
	 {
	  flag = 1;
	  for (j = 0; (j < m) && (flag == 1); j++)
	  {
	   if (i == basis[j])
	   {
		 flag = 0;
	   }
	  }
	  if (flag == 1)
	  {
		d[pos++] = i;
	  }
	 } // for



	} // set_d




	public static void set_Initial_A_aux()
	{

	  int i;
	  int j;
	  int k;


        for (i = 0; i < Initial_n; i++)
        {
            k =Initial_d[i];
          
            for (j = 0; j < m; j++)
            {
                Initial_A_aux.mat[j, i] = Initial_A.mat[j, k];
               
            } // for 2
        } // for1


    } 


	public static void simplex_algorithm()
	{
	   int i;
	   int j;
	   int k;
	   int optimal_flag;
	   int enter_id;
	   int exiting_id;
	   int itemp;
	   int basis_i;
	   double dtemp;
	   double min_value;
	   int count = 0;

	  optimal_flag = 0;

        Console.WriteLine(" m = " + m + ", n = " + n);

        Console.WriteLine("basis1");
        for (i = 0; i < m; i++)
            Console.Write(basis[i] + " ");
        Console.WriteLine();


        Console.WriteLine("A:");
        for (i = 0; i < m; i++)
        {
            for (j = 0; j < n; j++)
                Console.Write(A.mat[i, j] + " ");
            Console.WriteLine();
        } // for



        while (optimal_flag == 0)
        {

            count++;

            bublesort(basis, m);

            Console.WriteLine("basis2:");
            for (i = 0; i < m; i++)
                Console.Write(basis[i] + " ");
            Console.WriteLine();


            set_d();


            Console.WriteLine("d:");
            for (i = 0; i < n; i++)
                Console.Write(d[i] + " ");
            Console.WriteLine();


            set_A_aux();

            Console.WriteLine("A_aux (B, D):");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    Console.Write(A_aux.mat[i, j] + " ");
                Console.WriteLine();
            } // for

			B.copy_submatrix(A_aux, 0, m, 0, m);

            Console.WriteLine("B:");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                    Console.Write(B.mat[i, j] + " ");
                Console.WriteLine();
            } // for  

			BI.inv_gaussian(B, m);

            BI.erase_epsilons_matrix( m, m);

            Console.WriteLine("BI:");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < m; j++)
                    Console.Write(BI.mat[i, j] + " ");
                Console.WriteLine();
            } // for  


           // BIA_aux.mat= BI * A_aux;
            matrix_mult(BIA_aux, BI, A_aux,m, m, n);
			 BIA_aux.erase_epsilons_matrix( m, n);

            Console.WriteLine("BIA_aux (I, B-1*D):");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    Console.Write(BIA_aux.mat[i, j] + " ");
                Console.WriteLine();
            } // for  

            Console.WriteLine("A_aux (B,D):");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n; j++)
                    Console.Write(A_aux.mat[i, j] + " ");
                Console.WriteLine();
            } // for  


            Console.WriteLine("b:");
            for (i = 0; i < m; i++)
                Console.Write(b[i] + " ");
			Console.WriteLine();


            //BIb=BI* b;
			
            matrix_vector_mult(BIb, BI, b, m, m);
            erase_epsilons_vector(BIb ,m);
            Console.WriteLine("BIb:");
            for (i = 0; i < m; i++)
                Console.Write(BIb[i] + " ");
            Console.WriteLine();


            D.copy_submatrix( A_aux,0, m, m, n - m); 

            Console.WriteLine("D:");
            for (i = 0; i < m; i++)
            {
                for (j = 0; j < n - m; j++)
                    Console.Write(D.mat[i, j] + " ");
                Console.WriteLine();
            } // for  

            // END OF FOR DEBUG ONLY

            compute_cb_cd();

            Console.WriteLine("cb:");
            for (i = 0; i < m; i++)
                Console.Write(cb[i] + " ");
            Console.WriteLine();

            Console.WriteLine(" cd:");
            for (i = 0; i < (n - m); i++)
                Console.Write(cd[i] + " ");
            Console.WriteLine();


            // cbBI = cb * B-1

            vector_matrix_mult(cbBI, cb, BI, m, m);
            erase_epsilons_vector(cbBI, m);

            Console.WriteLine("cbBI:");
            for (i = 0; i < m; i++)
                Console.Write(cbBI[i] + " ");
            Console.WriteLine();


            // cbBID= cbBI* D;
            vector_matrix_mult(cbBID, cbBI, D,
                 m, n - m);
            erase_epsilons_vector(cbBID, n - m);


            Console.WriteLine("cbBID:");
            for (i = 0; i < (n - m); i++)
                Console.Write(cbBID[i] + " ");
            Console.WriteLine();



            vector_subtract(rd, cd, cbBID,
             n - m);
            erase_epsilons_vector(rd, n - m);

            Console.WriteLine("rd( cd - cbBID ):");
            for (i = 0; i < (n - m); i++)
                Console.Write(rd[i] + " ");
            Console.WriteLine();


            min_value = find_min_value(rd, n);
            if (min_value >= 0.0)
                optimal_flag = 1;
            else
            {
                enter_id = find_most_negative();
                exiting_id =BIA_aux.find_exiting_id( BIb, enter_id, n, m);
                Console.WriteLine("enter_id  = " + enter_id + ",  exiting_id = " + exiting_id +

                         " d[exiting_id] = " + d[exiting_id]);

                Console.WriteLine("pivot: enter_id = " + enter_id + ", exiting_id = " +
              d[exiting_id]);




                basis[exiting_id] = enter_id;
                Console.WriteLine("basis3:");
                for (i = 0; i < m; i++)
                    Console.Write(basis[i] + " ");
                Console.WriteLine();

            } // else 

        } // while
    } //simplex_algorithm

    public static void matrix_vector_mult(double[] c, Matrix A, double[] b, int n, int m)
    {
        int i, j, k;
        double sum;

        for (i = 0; i < n; i++)
        {
            sum = 0.0;
            for (k = 0; k < m; k++)
                sum = sum + A.mat[i, k] * b[k];
            c[i] = sum;
        } // for

    } // matrix_vector_mult

    public static void stage2_swap_colums(int i, int j)
	{
	  int k;
	  double temp;



	  for (k = 0; k < m; k++)
	  {
           // A_aux.
		 A_aux.mat[k,i] = A.mat[k,j];
		 A_aux.mat[k,j] = A.mat[k,i];
	  } // for

	} // swap_colums

	public static void swap_colums(double[][] A, int i, int j, int m, int n)
	{
	  int k;
	  double temp;

	  for (k = 0; k < m; k++)
	  {
		 temp = A[k][i];
		 A[k][i] = A[k][j];
		 A[k][j] = temp;
	  } // for


	} // swap_colums

	public static void swap_rows(Matrix W, int n, int m1, int m2)
	{
	  int i;
	  double temp;

	  for (i = 0; i <= 2 * n; i++)
	  {
		temp = W.mat[m1,i];
		W.mat[m1,i] = W.mat[m2,i];
		W.mat[m2,i] = temp;

	  } // for

	} // swap_rows


	public static void vector_matrix_mult(double[] c, double[] b, Matrix A, int n, int m)
	{
	  int i;
	  int j;
	  int k;
	  double sum;

	 for (i = 0; i < m; i++)
	 {
		   sum = 0.0;
		   for (k = 0; k < n; k++)
		   {
				   sum = sum + A.functorMethod(k,i) * b[k];
		   }
		 c[i] = sum;
	 } // for

	} // vector_matrix_mult


	public static void vector_subtract(double[] result_v, double[] v1, double[] v2, int n)
	{

	  int i;

	  for (i = 0; i < n; i++)
	  {
	   result_v[i] = v1[i] - v2[i];
	  }

	} // vector_subtract

	

}