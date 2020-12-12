package models;

import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;
import inf.v3d.view.Viewer;
import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.Node;
import fem_3D.Brick;
import fem_3D.BrickStructure;
import fem_3D.Visualizer3D;

public class SimpleBrick {
	public static BrickStructure create(){
		BrickStructure str = new BrickStructure(); 
		Node n1 = str.addNode(0.0, 0.0, 0.0);
		Node n2 = str.addNode(1.0, 0.0, 0.0);
		Node n3 = str.addNode(1.0, 1.0, 0.0);
		Node n4 = str.addNode(0.0, 1.0, 0.0);
		Node n5 = str.addNode(0.0, 0.0, 1.0);
		Node n6 = str.addNode(1.0, 0.0, 1.0);
		Node n7 = str.addNode(1.0, 1.0, 1.0);
		Node n8 = str.addNode(0.0, 1.0, 1.0);
		ShapeFunctions sf = new LinearSF("Brick");
		str.setSF(sf);
		
		double eModulus = 100;
		double nuy = 0.3333;
		Force f = new Force(2.0,0.0,0.0);
		
		Constraint c1 = new Constraint(false, false, false);
		Constraint c2 = new Constraint(true,  false, false);
		Constraint c3 = new Constraint(true,  true,  false);
		Constraint c4 = new Constraint(false, true,  false);
		Constraint c5 = new Constraint(false, false, true );
		Constraint c6 = new Constraint(true,  false, true );
		Constraint c7 = new Constraint(true,  true,  true );
		Constraint c8 = new Constraint(false, true,  true );
		
		n1.setConstraint(c1); n2.setConstraint(c2); n3.setConstraint(c3); n4.setConstraint(c4);
		n5.setConstraint(c5); n6.setConstraint(c6); n7.setConstraint(c7); n8.setConstraint(c8);
		
		n2.setForce(f);
		n3.setForce(f);
		n6.setForce(f);
		n7.setForce(f);
		
		Brick e1 = new Brick(eModulus, nuy, str.getLON());
		
		str.addElement(e1);
		str.assignStartIndex();
		return str;
	}
	
	public static void main (String [] args ) {
		
		BrickStructure brick = create();
		((Brick) brick.getLOE().get(0)).assignGaussPoints(2,2,2, brick.getSF());
		brick.solve();
		brick.printResults();
		brick.selectDisplacements();
		
		((Brick) brick.getElement(0)).computeStrain();
		((Brick) brick.getElement(0)).computeStress();
		
		Viewer viewer = new Viewer ();
		Visualizer3D viz = new Visualizer3D (brick,viewer);
		viz.setConstraintSymbolScale (0.08);
		viz.setForceSymbolScale (0.05);
		viz.setForceSymbolRadius (0.008);
		viz.setDisplacementScale (1);
		viz.drawElements ();
		viz.drawConstraints();
		viz.drawForces();
		viz.drawDisplacements();
		//viewer.setVisible(true);
	}
}