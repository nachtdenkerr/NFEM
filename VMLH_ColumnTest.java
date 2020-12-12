package models;

import java.util.ArrayList;
import java.util.List;

import fem_1D.Constraint;
import fem_1D.Force;
import fem_1D.GaussPoint;
import fem_1D.Node;
import plasticity.VMLH3D;
import shape_functions.LinearSF;
import shape_functions.ShapeFunctions;
import solver.NewtonRaphsonPla;
import fem_3D.Brick;
import fem_3D.BrickStructure;
import fem_abstract.Element;

public class VMLH_ColumnTest {
	@SuppressWarnings("unused")
	public static void main(String[] args) {
			
		BrickStructure struct = new BrickStructure();
		ShapeFunctions sf = new LinearSF("Brick");
		struct.setSF(sf);
	
		Force f = new Force(0, 1, 1);
	
		// create nodes
		Node n1 = struct.addNode(0, 0, 0);
		Node n2 = struct.addNode(1, 0, 0);
		Node n3 = struct.addNode(1, 1, 0);
		Node n4 = struct.addNode(0, 1, 0);
		Node n5 = struct.addNode(0, 0, 1);
		Node n6 = struct.addNode(1, 0, 1);
		Node n7 = struct.addNode(1, 1, 1);		
		Node n8 = struct.addNode(0, 1, 1);
		Node n9 = struct.addNode(0, 0, 2);
		Node n10= struct.addNode(1, 0, 2);
		Node n11= struct.addNode(1, 1, 2);		
		Node n12= struct.addNode(0, 1, 2);
		Node n13= struct.addNode(0, 0, 3);
		Node n14= struct.addNode(1, 0, 3);
		Node n15= struct.addNode(1, 1, 3);		
		Node n16= struct.addNode(0, 1, 3);
		
		Constraint c1 = new Constraint(false, false, false);
		
		// apply BCs
		n13.setForce(f);
		n14.setForce(f);
		n15.setForce(f);
		n16.setForce(f);
		
		n4.setConstraint(c1);
		n3.setConstraint(c1);
		n2.setConstraint(c1);
		n1.setConstraint(c1);
		
		// create elements
		double E = 1.0, nuy = 0.3;
		for (int i = 0; i < 3; i++){
			List<Node> ele = new ArrayList<Node>();
			for (int j = 0; j < 8; j++){
				ele.add(struct.getLON().get(j+4*i));
			}
			Element brick = new Brick (E, nuy, ele);
			((Brick) brick).assignGaussPoints(2, 2, 2, sf);
			struct.addElement(brick);
		}
		
		struct.assignStartIndex();
		struct.computeBMatrix(sf);
		
		for (Element b : struct.getLOE()){
			for (GaussPoint gp : b.getLOGP()){
				VMLH3D vm = new VMLH3D(gp, 5, 1.0/5.0, E,nuy);
				gp.setPlasticModel(vm);
			}
		}

		NewtonRaphsonPla returnMapping = new NewtonRaphsonPla(struct);
		returnMapping.solve(2, 0.0/10.0, 1.0/10.0);
	}
}