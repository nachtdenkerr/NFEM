package test;

import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.GaussPoint;
import fem_1D.Material;
import fem_1D.Node;
import fem_3D.Brick;
import fem_3D.BrickStructure;
import fem_abstract.Element;
import java.util.ArrayList;
import java.util.List;

import plasticity.PlasticModel;
import plasticity.VM3D;
import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;
import solver.NewtonRaphsonPla;

public class TwoBrickPatchTest_VM3D {
	public static void main (String args[]){
		
		BrickStructure bstr = new BrickStructure();
		ShapeFunctions sf = new LinearSF("Brick");
		bstr.setSF(sf);

		for (int i = 0; i <= 2; i++){
			bstr.addNode(0.0,0.0,1.0*i);
			bstr.addNode(1.0,0.0,1.0*i);
			bstr.addNode(1.0,1.0,1.0*i);
			bstr.addNode(0.0,1.0,1.0*i);
		}
		
		Constraint c1 = new Constraint(false, false, false);	
		for (int i = 0; i < 4; i++){
			bstr.getLON().get(i).setConstraint(c1);
		}
	
		double E = 250.0; double nuy = 0.3333;
		double yield_stress = 20.0;
		double hardening = 2.0*E/5.0;
		
		List<Node> ele1 = new ArrayList<Node>();
		List<Node> ele2 = new ArrayList<Node>();

		for (int i = 0; i < 8; i++){
			ele1.add(bstr.getLON().get(i));
			ele2.add(bstr.getLON().get(i+4));
		}
		
		Brick brick1 = new Brick (E, nuy, ele1);
		brick1.assignGaussPoints(2, 2, 2, sf);
		Brick brick2 = new Brick (E, nuy, ele2);
		brick2.assignGaussPoints(2, 2, 2, sf);

		bstr.addElement(brick1);
		bstr.addElement(brick2);
		
		bstr.assignStartIndex();
		
		Force f = new Force (0.0, 0.0, 20.0);
		bstr.listOfElements.get(1).getLON().get(2).setForce(f);
		bstr.listOfElements.get(1).getLON().get(3).setForce(f);
		
		System.out.println("************************************************************");
		System.out.println("*********************VON-MISES MATERIAL*********************");
		System.out.println("************************************************************");
		System.out.println();

		for (Element b : bstr.getLOE()){
			for (GaussPoint gp : b.getLOGP()){
				Material mat = new Material(2);
				gp.setMaterial(mat);
				PlasticModel vm = new VM3D(gp,yield_stress,hardening,E,nuy);
				gp.setPlasticModel(vm);
			}
		}
		bstr.computeBMatrix(sf);
		
		NewtonRaphsonPla returnMapping = new NewtonRaphsonPla(bstr);
		returnMapping.solve(6, 1.0/10.0, 1.0/10.0);
	}
}