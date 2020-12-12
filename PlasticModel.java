package plasticity;

import java.util.List;

import iceb.jnumerics.IVectorRO;
import iceb.plotutils.ui.PlotPanel;

public abstract class PlasticModel {
	@SuppressWarnings("rawtypes")
	public abstract List returnMapping(IVectorRO EPSN1, IVectorRO DEPS, IVectorRO EPSEN, double EPBARN, IVectorRO STRESN);
	
	public abstract void tangentUpdate (double EPBAR);

	public PlotPanel drawPlot() {
		// TODO Auto-generated method stub
		return null;
	}
}