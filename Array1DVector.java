package linalg;

import java.util.Arrays;

//import fem.Vector;
import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVector;
import iceb.jnumerics.IVectorRO;
import iceb.jnumerics.Vector3D;

@SuppressWarnings("serial")
public class Array1DVector implements IVectorRO,IVector{
	private double[] elements;

	public Array1DVector(double... x) {
		this.elements = x;
	}

	public Array1DVector(int n) {
		this.elements = new double[n];
	}

	public Array1DVector(Vector3D v) {
		for (int i = 0; i < v.getSize(); i++){
			this.elements = new double[v.getSize()];
			this.elements[i] = v.get(i);
		}
	}

	public IVectorRO crossProduct (IVectorRO direction){
		IVectorRO r = new Array1DVector(this.getSize());
		((Array1DVector) r).set(0, this.get(1)*direction.get(2) - this.get(2)*direction.get(1));
		((Array1DVector) r).set(1,-(this.get(0)*direction.get(2) - this.get(2)*direction.get(0)));
		((Array1DVector) r).set(2, this.get(0)*direction.get(1) - this.get(1)*direction.get(0));
		return r;
	}
	
	@Override
	public String toString() {
		return Arrays.toString(this.elements);
	}

	public void print(String l) {
		System.out.print(l + " = (");
		System.out.printf("%12.5e", get(0));
		for (int i = 1; i < this.elements.length; i++) {
			System.out.printf(", %12.5e", get(i));
		}
		System.out.println(")^T");
	}

	/*public Vector copy() {
		Vector c = new Vector(this.getSize());
		for (int i = 0; i < this.getSize(); i++) {
			c.set(i, this.get(i));
		}
		return c;
	}*/
	
	@Override
	public IVectorRO add(IVectorRO y) {
		IVectorRO r = new Array1DVector(this.getSize());
		
		for (int i = 0; i < this.getSize(); i++) {
			((Array1DVector) r).set(i, this.get(i) + y.get(i));
		}
		return r;
	}
	
	@Override
	public IVectorRO add(double alpha, IVectorRO y) {
		IVectorRO r = new Array1DVector(this.getSize());
		
		for (int i = 0; i < this.getSize(); i++) {
			((Array1DVector) r).set(i, this.get(i) + alpha * y.get(i));
		}
		return r;
	}

	@Override
	public void assignTo(double[] arg0) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public double dot(IVectorRO y) {
		double s = 0;

		for (int i = 0; i < this.getSize(); i++) {
			s += this.get(i) * y.get(i);
		}
		return s;
	}

	
	public IVectorRO multiplyMatrix(IMatrixRO a){
		int m = this.getSize();
		int n = a.getColumnCount();
		IVectorRO r = new Array1DVector(n);
		for (int i = 0; i < n; i++){
			double s = 0;
			for (int j = 0; j < m; j++){
				s = s + this.get(j)*a.get(j, i);
			}
			((Array1DVector) r).set(i, s);
		}
		return r;
	}
	
	public IMatrixRO dyadicProduct(IVectorRO y) {
		int m = this.getSize();
		int n = y.getSize();
		IMatrixRO r = new Array2DMatrix(m,n);
		for (int i = 0; i < m; i++){
			for (int j = 0; j < n; j++){
				((Array2DMatrix) r).set(i, j, this.get(i)*y.get(j));
			}
		}
		
		return r;
	}

	@Override
	public void dyadicProduct(double arg0, IVectorRO arg1, IMatrix arg2) {
		// TODO Auto-generated method stub
		
	}

	@Override
	public boolean equals(IVectorRO y) {
		int m = this.getSize();
		int n = y.getSize();
		if (m==n){
			for (int i = 0; i < m;){
				if (this.get(i) == y.get(i)) return true;
				else return false;
			}
		}
		return false;
	}

	@Override
	public double get(int idx) {
		if (idx >= 0 && idx < getSize()) {
			return this.elements[idx];
		} else {
			return Double.NaN;
		}
	}

	@Override
	public int getSize() {
		return this.elements.length;
	}

	@Override
	public IVectorRO multiply(double alpha) {
		IVectorRO r = new Array1DVector(this.getSize());
		int n = this.getSize();
		for (int i = 0; i < n; i++) {
			((Array1DVector) r).set(i, alpha * this.get(i));
		}
		return r;
	}

	@Override
	public IVector mutableCopy() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public double normInfinity() {
		double max = 0;
		for (int i = 0; i < this.getSize(); i++) {
			double v = Math.abs(this.get(i));
			if (max < v) {
				max = v;
			}
		}
		return max;
	}

	@Override
	public double normOne() {
		double s = 0;
		for (int i = 0; i < this.getSize(); i++) {
			s += Math.abs(this.get(i));
		}
		return s;
	}

	@Override
	public double normTwo() {
		return Math.sqrt(this.dot(this));
	}

	@Override
	public IVectorRO subtract(IVectorRO y) {
		IVectorRO r = new Array1DVector(this.getSize());
		
		for (int i = 0; i < this.getSize(); i++) {
			((Array1DVector) r).set(i, this.get(i) - y.get(i));
		}
		return r;
	}

	@Override
	public double[] toArray() {
		double[] r = new double[this.getSize()];
		for (int i = 0; i < r.length; i++){
			r[i] = this.get(i);
		}
		return r;
	}
	
	@Override
	public void add(int idx, double y) {
		for (int i = 0; i < this.getSize(); i++) {
			this.set(i, this.get(i));
		}
		this.set(idx, this.get(idx)+y);
	}

	@Override
	public void assignFrom(double[] v) {
		for (int i = 0; i < v.length; i++){
			this.set(i, v[i]);
		}
	}
	
	@Override
	public void assignFrom(IVectorRO v) {
		for (int i = 0; i < v.getSize() ; i++){
			this.set(i, v.get(i));
		}
	}

	@Override
	public IVector clone() {
		// TODO Auto-generated method stub
		return null;
	}

	@Override
	public void fill(double v) {
		Arrays.fill(this.elements, v);
	}

	@Override
	public void set(double... arg0) {
		for (int i = 0; i < arg0.length; i++){
			this.set(i, arg0[i]);
		}
	}

	@Override
	public void set(int idx, double v) {
		if (idx >= 0 && idx < getSize()) {
			this.elements[idx] = v;
		}
	}
}
