/******************************************************************************
 *  Compilation:  javac -d . -classpath C:/gurobi801/win64/lib/gurobi.jar;. TestSparseGurobiDCA_Qua.java
 *  Execution:    java -classpath C:/gurobi801/win64/lib/gurobi.jar;. TestSparseGurobiDCA_Qua
 *
 *  Ubuntu Compilation: javac -cp ".:$HOME/Documents/gurobi801/linux64/lib/gurobi.jar" TestSparseGurobiDCA_Qua.java
 *  Ubuntu Execution: java -cp ".:$HOME/Documents/gurobi801/linux64/lib/gurobi.jar" TestSparseGurobiDCA_Qua
 *  
 *  Please change the path above to the directory of gurobi
 *
 *  A test class to evaluate the numerical performance of DCA, and Gurobi for DRAP-NC 
 *  with quadratic objectives
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/


import java.util.*;
import java.io.*;
import gurobi.*;
import java.text.DecimalFormat;

public class TestSparseGurobiDCA_Qua{
	/**
     * An internal class to define the Objective functions
     * The name is QuadraticFunction, you can change it to any (convex) function 
     * by overwriting the abstract method getValue(double x)  
     * @param a, b, etc parameter of the function
     */
	//linear function
	private static class QuadraticFunction extends Function {
		double a;
		double b;
		public QuadraticFunction(double a, double b) {
			this.a = a;
			this.b = b;
		}
		public double getValue(double x) {
			return a * x * x + b * x;
		}
	}

	/**
     * main method
     * initial the test and record the test result in a text file
     * You can modify the name, the header of the text file as well as the test cases
     */
	public static void main(String[] args) throws FileNotFoundException{
		//Set output environment
		long currentTime = System.currentTimeMillis();
		PrintStream o = new PrintStream(new File("DCA----Gurobi" + " Quadratic + Sparse formulation" + Long.toString(currentTime) + "   Vb 10" + ".txt"));
		System.setOut(o);

		
		//test settings (dimension)
		//For gurobi, it is recommend to set size less than 10000.

		int[] testSize = new int[]{400, 800, 1600, 3200, 6400, 9600, 12800};
		//testSize = new int[]{400, 400<<1, 400<<2, 400<<3, 400<<4, 400<<5, 400<<6, 400<<7};
		//testSize = new int[]{12800, 25600, 51200};
		//int x = 51200;
		//testSize = new int[]{25600, 52100, 51200<<1, 51200<<2, 51200<<3, 51200<<4};
		//testSize = new int[]{400, 800, 1600, 3200, 6400, 9600, 12800, 25600, 52100, 51200<<1, 51200<<2, 51200<<3, 51200<<4, 51200<<5, 51200<<6};
		testSize = new int[]{50, 100, 200, 220, 240, 260, 280, 300, 320, 340, 360, 280, 300, 320, 340, 360, 380, 400};
		int rep = 10;
		int varBound = 10;

		System.out.println("Test: Variable bound " + varBound + ": objectives type: Quadratic Function + Sparse formulation");
		System.out.println("Dimension        DCA    Gurobi   Gap      NodeCount(Explored)    RootGap      HGap");

		for (int i = 0; i < testSize.length; i++) {
			for (int j = 0; j < rep; j++) {
				test(testSize[i], varBound);
			}
			System.out.println(" ");
		}
		
		return;
	}

	/**
     * test method
     * This is a static method that performs the tests.
     * You can modify the parameters to perform any tests 
     * The generation procedure ensures that every test instance is feasible
     * Default:
     * @param size dimension of the problem
     * @param varBound related to the upper bound of the box constraints: d_i < 1.3 * varBound
     */
    public static void test (int size, int varBound) {

    	//generate a feasible RAPNC instance
    	List<Function> obj = new ArrayList<Function>();		
		
		int dimension = size; 

		long[] lbVar = new long[dimension];
		long[] ubVar = new long[dimension];
		long[] lbNested = new long[dimension];
		long[] ubNested = new long[dimension];
		long B = 0;
		//for gurobi
		String[] varNames = new String[dimension];
		String[] varNamesY = new String[dimension];

    	try {
    		GRBEnv env = new GRBEnv("qp.log");    		
      		GRBModel model = new GRBModel(env);
      		model.set("OutputFlag", "1");   
      		model.set("TimeLimit", "1200.0");   		
      		//model.set(GRB.DoubleParam.FeasibilityTol, 0.000001);
      		double[] lbVarGurobi = new double[dimension];
      		double[] ubVarGurobi = new double[dimension];
      		char[] typeGurobi = new char[dimension];
			// Create variables
      		for (int i = 0; i < dimension; i++) {
				long b = (long) ((Math.random() + 0.3)* 10);
				long c = (long) (Math.random() * 10);
				if (b < c) {
					lbVar[i] = b;
					ubVar[i] = c;
				} else {
					lbVar[i] = c;
					ubVar[i] = b;
				}				
				lbVarGurobi[i] = (double) lbVar[i];
				ubVarGurobi[i] = (double) ubVar[i];		
				typeGurobi[i] = GRB.INTEGER;

				varNames[i] = "x" + Integer.toString(i);
		    }		
			
			GRBVar[] varlist = model.addVars(lbVarGurobi, ubVarGurobi, null, typeGurobi, varNames);
      		
      		// Set objective
			//linear objectives
      		GRBQuadExpr objGurobi = new GRBQuadExpr();
      		for (int i = 0; i < dimension; i++) {
      			double a = Math.random() * dimension;
      			double b = Math.random() * dimension;
      			
      			if (Math.random() < 0.5) {
      				b = -b;
      			}
      			QuadraticFunction fi = new QuadraticFunction(a, b);
				obj.add(fi);

				objGurobi.addTerm(b, varlist[i]);
				objGurobi.addTerm(a, varlist[i], varlist[i]);
      		}
            
      		model.setObjective(objGurobi);


      		//Nested constraints
      		double a = 0;
      		double b = 0;
      		double[] lbNestGurobi = new double[dimension];
      		double[] ubNestGurobi = new double[dimension];
      		//reset from 0 to dimension - 1
      		for (int i = 0; i < dimension - 1; i++) {
      			a += lbVar[i] + Math.random() * (ubVar[i] - lbVar[i]);
      			b += lbVar[i] + Math.random() * (ubVar[i] - lbVar[i]);
      			if (a > b) {
      				lbNested[i] = (long) b;
      				ubNested[i] = (long) a;
      			} else {
      				lbNested[i] = (long) a;
      				ubNested[i] = (long) b;
      			}

      			lbNestGurobi[i] = (double) lbNested[i];
      			ubNestGurobi[i] = (double) ubNested[i];
      			varNamesY[i] = "y" + Integer.toString(i); 
      			//GRBLinExpr exprNested = new GRBLinExpr();
      			//for (int j = 0; j <= i; j++) {
      			//	exprNested.addTerm(1.0, varlist[j]);
      			//}
      			//model.addConstr(exprNested, GRB.LESS_EQUAL, (double) ubNested[i], "u" + Integer.toString(i));
      			//model.addConstr(exprNested, GRB.GREATER_EQUAL, (double) lbNested[i], "l" + Integer.toString(i));
      		}

      		
            // Add constraint:\sum x_i = B
      		a += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      		b += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      		B = (long) Math.max(a, b);
      		// DCA
      		lbNested[dimension - 1] = B;
      		ubNested[dimension - 1] = B;
      		lbNestGurobi[dimension - 1] = (double) lbNested[dimension - 1];
      		ubNestGurobi[dimension - 1] = (double) ubNested[dimension - 1];

      		//GRBLinExpr expr = new GRBLinExpr();
      		//for (int i = 0; i < dimension; i++) {
      		//	expr.addTerm(1.0, varlist[i]);
      		//}
      		//model.addConstr(expr, GRB.EQUAL, B, "c0");

      		GRBVar[] varlistY = model.addVars(lbNestGurobi, ubNestGurobi, null, typeGurobi, varNamesY);
      		//nested constraint
      		for (int i = 0; i < dimension; i++) {
      			GRBLinExpr exprNested = new GRBLinExpr();
      			exprNested.addTerm(1.0, varlistY[i]);
      			exprNested.addTerm(-1.0, varlist[i]);
      			if (i > 0) {
      				exprNested.addTerm(-1.0, varlistY[i-1]);
      			}
      			model.addConstr(exprNested, GRB.EQUAL, 0, "NESTED" + Integer.toString(i));
      		}

      		//solve the problem by DCA algorithm
      		RAPNC test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
      		long startTime = System.currentTimeMillis();
			//ResultTypeRAPNC res = test.solveIntegerLinearDCA();
			ResultTypeRAPNC res = test.solveIntegerDCA();
			if (dimension < 10000) {
				for (int i = 0; i < 100; i++) {
					//res = test.solveIntegerLinearDCA();
					res = test.solveIntegerDCA();
				}
			}
			
			//time 
			long endTime = System.currentTimeMillis();
			long timeDCA = endTime - startTime;
			/*
			if (dimension < 10000) {
				timeDCA /= 20;
			}*/
    		//System.out.println("That took " + timeDCA + " milliseconds");

 			// Optimize Gurobi model
 			startTime = System.currentTimeMillis();
			model.optimize();
			endTime = System.currentTimeMillis();
			long timeGurobi = endTime - startTime;

			double MipGap = 0; //model.get(GRB.DoubleAttr.MIPGap);
			DecimalFormat df = new DecimalFormat("#.#####"); 
			String MipGapFormatted = df.format(MipGap); 
			int NodeCount = (int) model.get(GRB.DoubleAttr.NodeCount);
			double sum = 0;
      
			for (int i = 0; i < dimension; i++) {
				double val = varlist[i].get(GRB.DoubleAttr.X);
				if (Math.abs(val - res.sol[i]) >= 0.01) {
					System.out.println("No, the solution x[i] is different.");
					System.out.print(val);
					System.out.print(" ");
					System.out.println(res.sol[i]);
				}
				sum += obj.get(i).getValue(res.sol[i]);
			}

			if (Math.abs(sum - model.get(GRB.DoubleAttr.ObjVal)) > 0.001) {
				System.out.println(sum);
				System.out.println(model.get(GRB.DoubleAttr.ObjVal));
				System.out.println("No!!!!!!! The obj is not the same. But the values are close. Gurobi is not optimal");
			}

			// Root Gap 
			model.reset();
			model.set("NodeLimit", "0");      
			model.optimize();
			double RootGap = model.get(GRB.DoubleAttr.MIPGap);
			DecimalFormat rdf = new DecimalFormat("#.#####"); 
			String RootGapFormatted = rdf.format(RootGap);

			// Node 1 Gap 
			model.reset();
			model.set("NodeLimit", "1");      
			model.optimize();
			double HGap = model.get(GRB.DoubleAttr.MIPGap);
			DecimalFormat hdf = new DecimalFormat("#.#####"); 
			String HGapFormatted = hdf.format(HGap);

			System.out.println(
				String.format("%10s", dimension) 
				+ String.format("%10s", ((double) timeDCA) / 1000) 
				+ String.format("%10s", ((double) timeGurobi) / 1000) 
				+ String.format("%10s", MipGapFormatted) 
				+ String.format("%20s", NodeCount) 
				+ String.format("%10s", RootGapFormatted) 
				+ String.format("%10s", HGapFormatted) 
			);
      		// Dispose of model and environment 
            
      		model.dispose();
      		env.dispose();

    	} catch (GRBException e) {
      		System.out.println("Error code: " + e.getErrorCode() + ". " +
          	e.getMessage());
    	}
    }

}