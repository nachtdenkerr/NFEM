package solver;

import linalg.Array1DVector;
import iceb.jnumerics.*;

public class LUDecomposition {
	private static final double EPS = 1e-14;

	private IMatrix a;

	public LUDecomposition(IMatrixRO k_t) {
		initialize(k_t);
		decompose();
	}

	private void initialize(IMatrixRO k_t) {
		int n = k_t.getRowCount();

		this.a = new Array2DMatrix(n, n);
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				this.a.set(i, j, k_t.get(i, j));
			}
		}
	}

	private void decompose() {
		int n = this.a.getRowCount();

		for (int i = 0; i < n - 1; i++) {
			double aii = this.a.get(i, i);

			checkDiagonalElement(i, aii);
			for (int j = i + 1; j < n; j++) {
				double aji = this.a.get(j, i) / aii;
				this.a.set(j, i, aji);

				for (int k = i + 1; k < n; k++) {
					double aik = this.a.get(i, k);
					double ajk = this.a.get(j, k) - aik * aji;

					this.a.set(j, k, ajk);
				}
			}
		}
		checkDiagonalElement(n - 1, this.a.get(n - 1, n - 1));
	}

	private void checkDiagonalElement(int i, double aii) {
		if (Math.abs(aii) < EPS) {
			throw new IllegalArgumentException(
					"Zero diagonal element in row " + i);
		}
	}

	public IMatrix getL() {
		int n = this.a.getRowCount();
		IMatrix l = new Array2DMatrix(n, n);

		for (int j = 0; j < n; j++) {
			l.set(j, j, 1);
			for (int i = j + 1; i < n; i++) {
				l.set(i, j, this.a.get(i, j));
			}
		}
		return l;
	}

	public IMatrix getU() {
		int n = this.a.getRowCount();
		IMatrix u = new Array2DMatrix(n, n);

		for (int i = 0; i < n; i++) {
			for (int j = i; j < n; j++) {
				u.set(i, j, this.a.get(i, j));
			}
		}
		return u;
	}

	public IVectorRO solveFor(IVectorRO res) {
		IVectorRO y = computeY(res);
		IVectorRO x = computeX(y);

		return x;
	}

	public IVectorRO computeX(IVectorRO y) {
		int n = this.a.getRowCount();
		IVectorRO x = new Array1DVector(n);

		((Array1DVector) x).set(n - 1, y.get(n - 1) / this.a.get(n - 1, n - 1));
		for (int i = n - 2; i >= 0; i--) {
			double s = 0;

			for (int j = i + 1; j < n; j++) {
				s += this.a.get(i, j) * x.get(j);
			}
			((Array1DVector) x).set(i, (y.get(i) - s) / this.a.get(i, i));
		}
		return x;
	}

	public IVectorRO computeY(IVectorRO res) {
		int n = this.a.getRowCount();
		IVectorRO y = new Array1DVector(n);

		((Array1DVector) y).set(0, res.get(0));
		for (int i = 1; i < n; i++) {
			double s = 0;

			for (int j = 0; j < i; j++) {
				s += this.a.get(i, j) * y.get(j);
			}
			((Array1DVector) y).set(i, res.get(i) - s);
		}
		return y;
	}
}