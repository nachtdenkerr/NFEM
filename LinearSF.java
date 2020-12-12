package shape_functions;

import iceb.jnumerics.IMatrix;
import iceb.jnumerics.IMatrixRO;

import java.util.ArrayList;
import java.util.List;

import fem_1D.GaussPoint;
import linalg.Array2DMatrix;

public class LinearSF extends ShapeFunctions{
	String element_name;
	double eta1,eta2,eta3;
	public List<Double> listOfSFs = new ArrayList<>();
	public IMatrix d1sf = new Array2DMatrix (1,2);
	public IMatrix d2sf = new Array2DMatrix (2,4);
	public IMatrix dcst = new Array2DMatrix (2,3);
	public IMatrix d3sf = new Array2DMatrix (3,8);

	
	public LinearSF(String e){
		this.element_name = e;
	}
	
	public void setEta(GaussPoint gp){
		this.eta1 = gp.getPosition().get(0);
		this.eta2 = gp.getPosition().get(1);
		this.eta3 = gp.getPosition().get(2);
	}
	
	public List<Double> getSF(GaussPoint g){
		double a1 = 0.5*(1-eta1); double a2 = 0.5*(1+eta1);
		double b1 = 0.5*(1-eta2); double b2 = 0.5*(1+eta2);
		double c1 = 0.5*(1-eta3); double c2 = 0.5*(1+eta3);
		setEta(g);
		if (element_name == "Truss"){
			double sf_n1 = 0.5*(1-eta1);
			double sf_n2 = 0.5*(1+eta1);
			listOfSFs.add(sf_n1); listOfSFs.add(sf_n2);
		}
		
		if (element_name == "Quad"){
			double sf_n1 = a1*b1;
			double sf_n2 = a2*b1;
			double sf_n3 = a2*b2;
			double sf_n4 = a1*b2;
			listOfSFs.add(sf_n1); listOfSFs.add(sf_n2); listOfSFs.add(sf_n3); listOfSFs.add(sf_n4);
		}
		
		if (element_name == "CST"){
			double sf_n1 = eta1;
			double sf_n2 = eta2;
			double sf_n3 = 1 - eta1 - eta2;
			listOfSFs.add(sf_n1); listOfSFs.add(sf_n2); listOfSFs.add(sf_n3);
		}
		
		if (element_name == "Brick"){
			double sf_n1 = a1*b1*c1;
			double sf_n2 = a2*b1*c1;
			double sf_n3 = a2*b2*c1;
			double sf_n4 = a1*b2*c1;
			double sf_n5 = a1*b1*c2;
			double sf_n6 = a2*b1*c2;
			double sf_n7 = a2*b2*c2;
			double sf_n8 = a1*b2*c2;
			listOfSFs.add(sf_n1); listOfSFs.add(sf_n2); listOfSFs.add(sf_n3); listOfSFs.add(sf_n4);
			listOfSFs.add(sf_n5); listOfSFs.add(sf_n6); listOfSFs.add(sf_n7); listOfSFs.add(sf_n8);
		}
		
		return listOfSFs;
	}
	
	public IMatrixRO getdSF(GaussPoint g){
		setEta(g);
		IMatrixRO dsf = null;
		if (element_name == "Truss"){
			dsf = new Array2DMatrix (1,2);
			d1sf.set(0, 0, -0.5);
			d1sf.set(0, 1, 0.5);
			dsf = d1sf;
		}
		
		if (element_name == "Quad"){
			dsf = new Array2DMatrix (2,4);
			d2sf.set(0, 0, -0.25*(1-eta2));
			d2sf.set(0, 1, 0.25*(1-eta2));
			d2sf.set(0, 2, 0.25*(1+eta2));
			d2sf.set(0, 3, -0.25*(1+eta2));
			d2sf.set(1, 0, -0.25*(1-eta1));
			d2sf.set(1, 1, -0.25*(1+eta1));
			d2sf.set(1, 2, 0.25*(1+eta1));
			d2sf.set(1, 3, 0.25*(1-eta1));
			dsf = d2sf;
		}
		
		if (element_name == "CST"){
			dsf = new Array2DMatrix (2,3);
			dcst.set(0, 0, 1);
			dcst.set(0, 1, 0);
			dcst.set(0, 2, -1);
			dcst.set(1, 0, 0);
			dcst.set(1, 1, 1);
			dcst.set(1, 2, -1);
			dsf = dcst;
		}
		
		double a1 = 0.5*(1-eta1); double a2 = 0.5*(1+eta1);
		double b1 = 0.5*(1-eta2); double b2 = 0.5*(1+eta2);
		double c1 = 0.5*(1-eta3); double c2 = 0.5*(1+eta3);
		
		if (element_name == "Brick"){
			dsf = new Array2DMatrix (3,8);
			d3sf.set(0, 0,-0.5*b1*c1);
			d3sf.set(0, 1, 0.5*b1*c1);
			d3sf.set(0, 2, 0.5*b2*c1);
			d3sf.set(0, 3,-0.5*b2*c1);
			d3sf.set(0, 4,-0.5*b1*c2);
			d3sf.set(0, 5, 0.5*b1*c2);
			d3sf.set(0, 6, 0.5*b2*c2);
			d3sf.set(0, 7,-0.5*b2*c2);
			
			d3sf.set(1, 0,-0.5*a1*c1);
			d3sf.set(1, 1,-0.5*a2*c1);
			d3sf.set(1, 2, 0.5*a2*c1);
			d3sf.set(1, 3, 0.5*a1*c1);
			d3sf.set(1, 4,-0.5*a1*c2);
			d3sf.set(1, 5,-0.5*a2*c2);
			d3sf.set(1, 6, 0.5*a2*c2);
			d3sf.set(1, 7, 0.5*a1*c2);
			
			d3sf.set(2, 0,-0.5*b1*a1);
			d3sf.set(2, 1,-0.5*b1*a2);
			d3sf.set(2, 2,-0.5*b2*a2);
			d3sf.set(2, 3,-0.5*b2*a1);
			d3sf.set(2, 4, 0.5*b1*a1);
			d3sf.set(2, 5, 0.5*b1*a2);
			d3sf.set(2, 6, 0.5*b2*a2);
			d3sf.set(2, 7, 0.5*b2*a1);
			dsf = d3sf;
		}
		
		return dsf;
	}
}
