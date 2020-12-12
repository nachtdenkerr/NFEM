package fem_1D;

import iceb.jnumerics.IMatrixRO;
import iceb.jnumerics.IVectorRO;

public class Material {
	private IVectorRO epsE, sigma;
	private IMatrixRO tangent;
	private double epsP;
	private int hardeningLaw;
	
	public Material (int i) {
		if (i != 1 && i!= 2) {
			System.out.println("Wrong hardening law");
		}
		else this.hardeningLaw = i;

	}

	public void setEpsElastic(IVectorRO epsE){
		this.epsE = epsE;
	}
	
	public IVectorRO getEpsElastic(){
		return this.epsE;
	}
	
	public void setEpsPlastic(double epsP){
		this.epsP = epsP;
	}
	
	public double getEpsPlastic(){
		return this.epsP;
	}
	
	public void setSigma(IVectorRO sigma){
		this.sigma = sigma;
	}
	
	public IVectorRO getSigma(){
		return this.sigma;
	}
	
	public void setTangent(IMatrixRO tangent){
		this.tangent = tangent;
	}
	
	public IMatrixRO getTangent(){
		return this.tangent;
	}
	
	public double PLFUN(double sigma_y0, double H, double EPBAR) {
		double SIGMAY = -100;
		if (this.hardeningLaw == 1) {
			SIGMAY = sigma_y0 + H*(EPBAR);
		}
		else {
			SIGMAY = sigma_y0 + H*EPBAR*EPBAR;
//			SIGMAY = sigma_y0 + H*Math.sqrt(1.0 + EPBAR);
//			SIGMAY = sigma_y0 + H*(1.0 - Math.exp(-20.0*EPBAR));
		}
		return SIGMAY;
	}
	
	public double DPLFUN(double sigma_y0, double H, double EPBAR) {
		double HSLOPE = -100;
		if (this.hardeningLaw == 1) {
			HSLOPE = H;
		}
		else {
			HSLOPE = H*2.0*EPBAR;
//			HSLOPE = H*0.5*Math.pow(1.0 + EPBAR, -0.5);
//			HSLOPE = H*(-20.0)*(-Math.exp(-2.0*EPBAR));
		}
		return HSLOPE;
	}
}
