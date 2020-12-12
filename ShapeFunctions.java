package shape_functions;

import iceb.jnumerics.IMatrixRO;

import java.util.List;

import fem_1D.GaussPoint;

public abstract class ShapeFunctions {
	
	public abstract void setEta(GaussPoint gp);
	
	public abstract List<Double> getSF(GaussPoint g);
	
	public abstract IMatrixRO getdSF(GaussPoint g);
}