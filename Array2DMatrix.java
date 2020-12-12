package linalg;

import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.SolveFailedException;

@SuppressWarnings("serial")
public class Array2DMatrix implements IMatrixRO,IMatrix{
	private double[][] elements;

	public Array2DMatrix(int m, int n) {
		this.elements = new double[m][n];
	}

	public Array2DMatrix(int m, int n, double... values) {
		this(m, n);
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				this.set(i, j, values[i * n + j]);
			}
		}
	}

	public void set(int i, int j, double v) {
		this.elements[i][j] = v;
	}

	public void setRow(int idx, double... values) {
		for (int i = 0; i < this.getColumnCount(); i++) {
			this.set(idx, i, values[i]);
		}
	}

	public double get(int i, int j) {
		return this.elements[i][j];
	}

	public void print() {
		int m = getRowCount();
		//System.out.println(l + " = ");
		for (int i = 0; i < m; i++) {
			System.out.print("|");
			for (int j = 0; j < getColumnCount(); j++) {
				System.out.printf(" %12.5e", get(i, j));
			}
			System.out.println(" |");
		}
	}

	public Array2DMatrix transpose() {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		Array2DMatrix r = new Array2DMatrix(n, m);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				r.set(j, i, this.get(i, j));
			}
		}
		return r;
	}
	
	public double determinant(){
		int m = this.getRowCount();
		double det = 0;

        // Trivial 1x1 matrix
        if (m == 1) {
        	det = this.get(0, 0);
        	return det;
        }
        // Trivial 2x2 matrix
        else if (m == 2) {
        	det = this.get(0, 0)*this.get(1, 1) - this.get(1, 0)*this.get(0, 1);
        	return det;
        }
        
    	for (int i =0; i< m; i++){
    		IMatrixRO temp = new Array2DMatrix(m-1,m-1);
			for (int j = 1; j < m; j++) {
				for (int k = 0; k < m; k++) {
					if (k < i) {
						((Array2DMatrix) temp).set(j-1,k,this.get(j,k));
					} else if (k > i) {
						((Array2DMatrix) temp).set(j-1, k-1, this.get(j, k));
					}
				}
			}

			det = det + this.get(0,i)*Math.pow (-1,(double) i)*((Array2DMatrix) temp).determinant();
    	}
		return det;
	}

	@Override
	public IMatrixRO add(IMatrixRO b) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		IMatrixRO r = new Array2DMatrix(m,n);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				((Array2DMatrix) r).set(i, j, this.get(i, j) + b.get(i, j));
			}
		}
		return r;
	}

	@Override
	public IMatrixRO add(double alpha, IMatrixRO b) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		IMatrixRO r = new Array2DMatrix(m,n);
		
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				((Array2DMatrix) r).set(i, j, this.get(i, j) + alpha* b.get(i, j));
			}
		}
		return r;
	}

	@Override
	public void assignTo(double[][] arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public void assignTo(IMatrix arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean equals(IMatrixRO b) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				if (this.get(i, j) == b.get(i, j)){
					return true;
				}
			}
		}
		return false;
	}

	@Override
	public int getColumnCount() {
		int n = this.elements[0].length;
		return n;
	}

	@Override
	public int getRowCount() {
		int m = this.elements.length;
		return m;
	}
	
	public IMatrixRO getMatrix(int row_index, int col_index, int row_num, int col_num){
		int m = this.getRowCount();
		int n = this.getColumnCount();
		if (row_index + row_num > m || col_index + col_num > n){
			System.out.println("Unappropriate parameters");
			return null;
		}
		else {
			IMatrixRO r = new Array2DMatrix(row_num, col_num);
			return r;
		}
		
	}

	@Override
	public IMatrixRO invert() throws SolveFailedException {
		//Create a matrix to receive inverse matrix data
		IMatrixRO reMatrix = new Array2DMatrix(this.getRowCount(),this.getColumnCount()); 
		 //Get the value of the determinant of the original matrix
		 double value = this.determinant(); 
		 //Judge whether the value of matrix determinant is zero
		 if(Math.abs(value) <= 10e-6) {
		     System.out.println("The matrix is irreversible!");
		     return null;
		 }
		 //Primitive matrix mat Assignment divided by the value of the original determinant value Inverse matrix
		 for (int i = 0; i < reMatrix.getRowCount(); i++) {
		     for (int j = 0; j < reMatrix.getColumnCount(); j++) {
		         ((Array2DMatrix) reMatrix).set(i, j, this.getWithMatrix().get(i, j) / value);
		     }
		 }
		 return reMatrix;
	}
	
	public IMatrixRO getComplementMinor(IMatrixRO mat, int i, int j) {
		//Create a new matrix to receive the remainder expression, and delete the value of this column in this line
		IMatrixRO m = new Array2DMatrix(mat.getRowCount()-1,mat.getColumnCount()-1); 
		//To traverse a new matrix m Variables
		 int row =0 ,col=0;
		/*
		* Traversing the data of the original matrix, j2 represents row, k represents column
		*/
		for (int j2 = 0; j2 < mat.getRowCount(); j2++) {
		    //In the first place i Row division data omitted
			if(j2 == i) continue; 
	        for (int k = 0; k < mat.getColumnCount(); k++) {
	             //In the first place j Column data omitted
	             if(k == j) continue;
	             //assignment
	             ((IMatrix) m).set(row, col,mat.get(j2, k));
	             //Variables traversing a new matrix
	             col++;
	             if(col >= m.getColumnCount() ) {
	                 col = 0;
	                 row++;
	             }
	        }
	    }
	    return m;
	}
	
    public IMatrixRO getWithMatrix() {
         //Create a matrix to store the values of adjoint matrix
    	 IMatrixRO withMatrix = new Array2DMatrix(this.getRowCount(),this.getColumnCount());
         //ergodic withMatrix Store corresponding mat Value
         for (int i = 0; i < withMatrix.getRowCount(); i++) {
             for (int j = 0; j < withMatrix.getColumnCount(); j++) {
            	 IMatrixRO minor = getComplementMinor(this, j, i);
                 double temp = Math.pow(-1, i+j) *  ((Array2DMatrix) minor).determinant();
                 if(Math.abs(temp) <= 10e-6) temp = 0;
                 ((IMatrix) withMatrix).set(i, j,temp);
             }
         }
         //Return result
         return withMatrix;    
    }

	@Override
	public IMatrixRO multiply(double alpha) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		IMatrixRO r = new Array2DMatrix(m,n);
		
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				((Array2DMatrix) r).set(i, j, alpha*this.get(i, j));
			}
		}
		return r;
	}

	@Override
	public IMatrixRO multiply(IMatrixRO b) {
		int m = this.getRowCount();
		IMatrixRO r = new Array2DMatrix(m, b.getColumnCount());

		for (int i = 0; i < this.getRowCount(); i++) {
			for (int j = 0; j < ((Array2DMatrix) b).getColumnCount(); j++) {
				double s = 0;

				for (int k = 0; k < this.getColumnCount(); k++) {
					s += this.get(i, k) * b.get(k, j);
				}
				((Array2DMatrix) r).set(i, j, s);
			}
		}
		return r;
	}

	@Override
	public IVectorRO multiply(IVectorRO x) {
		IVectorRO r = new Array1DVector(this.getRowCount());

		for (int i = 0; i < this.getRowCount(); i++) {
			double s = 0;

			for (int j = 0; j < this.getColumnCount(); j++) {
				s = s + this.get(i, j) * x.get(j);
			}
			
			((Array1DVector) r).set(i, s);
		}
		return r;
	}

	@Override
	public IMatrixRO multiply(double alpha, IMatrixRO b) {
		Array2DMatrix r = new Array2DMatrix(this.getRowCount(), ((Array2DMatrix) b).getColumnCount());

		for (int i = 0; i < this.getRowCount(); i++) {
			for (int j = 0; j < ((Array2DMatrix) b).getColumnCount(); j++) {
				double s = 0;

				for (int k = 0; k < this.getColumnCount(); k++) {
					s += this.get(i, k) * b.get(k, j);
				}
				r.set(i, j, alpha * s);
			}
		}
		return r;
	}

	@Override
	public IVectorRO multiply(double alpha, IVectorRO x) {
		IVectorRO r = new Array1DVector(this.getRowCount());

		for (int i = 0; i < this.getRowCount(); i++) {
			double s = 0;

			for (int j = 0; j < this.getColumnCount(); j++) {
				s += this.get(i, j) * x.get(j);
			}
			((Array1DVector) r).set(i, alpha * s);
		}
		return r;
	}
	
	public IMatrixRO inverse2x2(){
		double det = 0;
		IMatrixRO inv = new Array2DMatrix(2,2);
		if (this.getColumnCount() == 2 && this.getColumnCount() == 2){
			det = this.get(0, 0)*this.get(1, 1) - this.get(0, 1)*this.get(1, 0);
			IMatrixRO adj = new Array2DMatrix(2,2);
			((Array2DMatrix) adj).set(0, 0, this.get(1,1));
			((Array2DMatrix) adj).set(1, 1, this.get(0, 0));
			((Array2DMatrix) adj).set(0, 1, -this.get(0,1));
			((Array2DMatrix) adj).set(1, 0, -this.get(1, 0));
			inv = adj.multiply(1/det);
		}
		return inv;
	}
	
	public IMatrixRO inverse3x3(){
		double det = this.determinant();
		IMatrix inv = new Array2DMatrix(3,3);
		double invdet = 1/det;
		inv.set(0, 0, (this.get(1, 1)*this.get(2, 2) - this.get(2, 1)*this.get(1, 2))*invdet);
		inv.set(0, 1, (this.get(0, 2)*this.get(2, 1) - this.get(0, 1)*this.get(2, 2))*invdet);
		inv.set(0, 2, (this.get(0, 1)*this.get(1, 2) - this.get(0, 2)*this.get(1, 1))*invdet);
		inv.set(1, 0, (this.get(1, 2)*this.get(2, 0) - this.get(1, 0)*this.get(2, 2))*invdet);
		inv.set(1, 1, (this.get(0, 0)*this.get(2, 2) - this.get(0, 2)*this.get(2, 0))*invdet);
		inv.set(1, 2, (this.get(1, 0)*this.get(0, 2) - this.get(0, 0)*this.get(1, 2))*invdet);
		inv.set(2, 0, (this.get(1, 0)*this.get(2, 1) - this.get(2, 0)*this.get(1, 1))*invdet);
		inv.set(2, 1, (this.get(2, 0)*this.get(0, 1) - this.get(0, 0)*this.get(2, 1))*invdet);
		inv.set(2, 2, (this.get(0, 0)*this.get(1, 1) - this.get(1, 0)*this.get(0, 1))*invdet);
		return inv;
	}

	@Override
	public IMatrix mutableCopy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double normInfinity() {
		double mn = Double.NEGATIVE_INFINITY;

		for (int i = 0; i < this.getRowCount(); i++) {
			for (int j = 0; j < this.getColumnCount(); j++) {
				double v = Math.abs(this.get(i, j));

				if (v > mn) {
					mn = v;
				}
			}
		}
		return mn;
	}

	@Override
	public double normOne() {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public double quadraticForm(IVectorRO arg0, IVectorRO arg1) {
		// TODO Auto-generated method stub
		return 0;
	}

	@Override
	public IVectorRO solve(IVectorRO arg0) throws SolveFailedException {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public IMatrixRO subtract(IMatrixRO b) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		Array2DMatrix r = new Array2DMatrix(m, n);

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				r.set(i, j, this.get(i, j) - b.get(i, j));
			}
		}
		return r;
	}

	public void add(int i, int j, double a) {
		this.elements[i][j] = this.elements[i][j] + a;
	}

	public void addColumn(int col_index, IVectorRO v) {
		int m = this.getRowCount();
		//int n = this.getColumnCount();
		for (int i = 0; i < m; i++) {
			this.set(col_index, i, this.get(col_index, i)+v.get(i));
		}
	}

	public void addColumn(int arg0, int arg1, IVectorRO arg2) {
		// TODO Auto-generated method stub
		
	}

	public void addMatrix(int row_index, int col_index, IMatrixRO b) {
		// TODO Auto-generated method stub
		int m = b.getRowCount();
		int n = b.getColumnCount();
		for (int i = row_index; i < row_index + m; i++) {
			for (int j = col_index; j < col_index + n; j++){
				this.set(i, j, this.get(i, j) + b.get(i-row_index, j - col_index));
			}
		}
	}

	public void addRow(int row_index, IVectorRO v) {
		//int m = this.getRowCount();
		int n = this.getColumnCount();
		for (int j = 0; j < n; j++) {
			this.set(row_index, j, this.get(row_index, j)+v.get(j));
		}
	}

	public void addRow(int a, int b, IVectorRO v) {
		
	}

	public void assignFrom(IMatrixRO arg0) {
		// TODO Auto-generated method stub
		
	}

	public boolean canWrite(int arg0, int arg1) {
		// TODO Auto-generated method stub
		return false;
	}

	@Override
	public IMatrix clone() {
		// TODO Auto-generated method stub
		return null;
	}

	public void fill(double a) {
		int m = this.getRowCount();
		int n = this.getColumnCount();

		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				this.elements[i][j] = a;
			}
		}
	}

	public void fillColumn(int col_index, double alpha) {
		int m = this.getRowCount();
		//int n = this.getColumnCount();
		for (int i = 0; i < m; i++) {
			this.set(i,col_index,alpha);
		}
	}

	public void fillDiagonal(double value) {
		int m = this.getRowCount();
		int n = this.getColumnCount();
		if (m == n){
			for (int i = 0; i < m; i++) {
				this.elements[i][i] = value;
			}
		}
		
	}

	public void fillRow(int row_index, double alpha) {
		//int m = this.getRowCount();
		int n = this.getColumnCount();
		for (int j = 0; j < n; j++) {
			this.set(row_index,j,alpha);
		}
	}

	@Override
	public IVector getColumn(int col_index) {
		int m = this.getRowCount();
		IVector result = new Array1DVector(m);
		for (int i = 0; i< m; i++){
			((Array1DVector) result).set(i, this.get(i, col_index));
		}
		return result;
	}

	@Override
	public IVector getRow(int row_index) {
		// TODO Auto-generated method stub
		int n = this.getColumnCount();
		IVector result = new Array1DVector(n);
		for (int j = 0; j< n; j++){
			((Array1DVector) result).set(j, this.get(row_index, j));
		}
		return result;
	}

	public void setColumn(int col_index, IVectorRO vector) {
		int m = this.getRowCount();
		if (vector.getSize() <= m){
			for (int i = 0; i< vector.getSize(); i++){
				this.set(i, col_index, vector.get(i));
			}
		}
	}

	public void setColumn(int row_index, int col_index, IVectorRO v) {
		int l = v.getSize();
		for (int i = 0; i < l; i++) {
			this.set(i + row_index, col_index, v.get(i));
		}
	}

	public void setMatrix(int row_index, int col_index, IMatrixRO a) {
		int m = a.getRowCount();
		int n = a.getColumnCount();
		for (int i = 0; i < m ; i++){
			for (int j = 0; j < n; j++){
				this.set(row_index + i, col_index + j, a.get(i , j));
			}
		}
	}

	public void setRow(int row_index, IVectorRO vector) {
		int n = this.getColumnCount();
		if (vector.getSize() <= n){
			for (int j = 0; j< vector.getSize(); j++){
				this.set(row_index, j, vector.get(j));
			}
		}
	}

	public void setRow(int row_index, int col_index, IVectorRO v) {
		int l = v.getSize();
		for (int j = 0; j < l; j++) {
			this.set(row_index, j + col_index, v.get(j));
		}
	}
	
	public IMatrixRO getIMatrix(int i){
		IMatrixRO I = new Array2DMatrix(i,i);
		for (int j = 0; j < i; j++){
			((Array2DMatrix) I).set(j,j,1);
		}
		return I;
	}

	@Override
	public void fillSubmatrix(int arg0, int arg1, int arg2, int arg3,
			double arg4) {
		// TODO Auto-generated method stub
		
	}
}
