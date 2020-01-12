/******************************************************************************
 *  Compilation:  javac -cp ./ TestFlow1Con.java
 *  Execution:    java -cp ./ TestFlow1Con
 *  
 *  A test class to evaluate the numerical performance of DCA, MDA, and FlowS
 *  @author Zeyang Wu @ University of Minnesota
 *
 ******************************************************************************/
import java.util.*;
import java.io.*;

public class TestFlow1Con {

	/**
     * An internal class to define the Objective functions
     * The name is ObjectiveFunction, you can change it to any (convex) function 
     * by overwriting the abstract method getValue(double x)  
     * @param a, b, etc. parameter of the function
     */
	private static class ObjectiveFunction extends Function {
		double a;
		double b;
		public ObjectiveFunction(double a, double b) {
			this.a = a;
			this.b = b;
		}
		public double getValue(double x) {
			//return b * x;
			//return a * x * x + b * x; //quadratic function
			//return 10*b + a / x; //Crash function
			return a * b * b / x / x / x; // Fuel function
			//return x * x * x * x / 4 +  b * x ; // F function
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

		//Please change the ObjectiveFunction above accordingly
		//String obj_function_type = " Linear ";
		//String obj_function_type = " Quadratic ";
		//String obj_function_type = " Con_Crash ";
		String obj_function_type = " Con_Fuel ";
		//String obj_function_type = " Con_F ";

		int rep = 10;
		int varBound = 100;

		PrintStream o = new PrintStream(new File(
			"DCA---MDA numerical experiment" 
			+ obj_function_type 
			+ Long.toString(currentTime) 
			+ "_Vb_ " 
			+ String.valueOf(varBound) + ".txt")
		);
		System.setOut(o);

		
		//test settings (dimension)
		int[] testSize = new int[]{800, 1600, 3200, 6400, 12800, 25600, 51200, 51200<<1};
		testSize = new int[]{800, 1600, 3200, 6400, 12800, 25600, 51200, 51200<<1, 51200<<2, 51200<<3, 51200<<4,  51200<<5, 51200<<6, 51200<<7};
		//testSize = new int[]{51200<<7, 51200<<8};
		//testSize = new int[]{52100<<2, 51200<<3, 51200<<4,  51200<<5, 51200<<6, 51200<<7, 51200<<8};
		//int x = 51200;
		//testSize = new int[]{5};
		

		System.out.println("Test: Variable bound " + varBound + ":General Convex Function");
		System.out.println("Dimension        DCA       FastMDA	      number_subproblem_DCA      number_subproblem_MDA");

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
    	List<Function> obj = new ArrayList<Function>();
		int dimension = size; 
		long[] lbVar = new long[dimension];
		long[] ubVar = new long[dimension];
		long[] lbNested = new long[dimension];
		long[] ubNested = new long[dimension];
		long B = (long) (0.75 * dimension * dimension);
		String[] varNames = new String[dimension];


		// Create box constraints on variables
      	for (int i = 0; i < dimension; i++) {
			long b = (long) ((Math.random() + 0.3)* varBound);
			long c = (long) (Math.random() * varBound);
			if (b < c) {
				lbVar[i] = b;
				ubVar[i] = c;
			} else {
				lbVar[i] = c;
				ubVar[i] = b;
			}				
			lbVar[i] = 0;
		}		
			
      		
      	// Set objective      		
      	for (int i = 0; i < dimension; i++) {
      		double a = Math.random();
      		double b = Math.random();
      		if (Math.random() < 0.5) {
      			b = b;
      		}
      		ObjectiveFunction fi = new ObjectiveFunction(a, b);
			obj.add(fi);
      	}
            

      	//Nested constraints
      	double a = 0;
      	double b = 0;
      	//reset from 1 to dimension
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
      	}
            
        // Add constraint:\sum x_i = B
      	a += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      	b += lbVar[dimension - 1] + Math.random() * (ubVar[dimension - 1] - lbVar[dimension - 1]);
      	B = (long) Math.max(a, b);
      	// Greedy Algorithm
      	lbNested[dimension - 1] = B;
      	ubNested[dimension - 1] = B;
      	
      	/*
      	System.out.println("ubNested");        
        System.out.println(Arrays.toString(ubNested));
        System.out.println("lbNested");
        System.out.println(Arrays.toString(lbNested));	
      	*/
      	
        /********************************************************************************
        **
        **DCA: solve the problem by DCA and record the time
        **
        ********************************************************************************/
      	//
      	RAPNC test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
      	long startTime = System.currentTimeMillis();
		ResultTypeRAPNC res = test.solveIntegerDCA();
		//System.out.println("DCA Done!");
		//ResultTypeRAPNC res = test.solveIntegerLinearDCA();
		//solve the problem by DCA 100 times
		for (int i = 0; i < 0; i++) {
			test = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
			//res = test.solveIntegerDCALinear();
			res = test.solveIntegerLinearDCA();
		}	

		//time 
		long endTime = System.currentTimeMillis();
		long timeDCA = endTime - startTime;
		//number of subproblems
		long number_subproblem_DCA = test.number_subproblem;
		/*
    	System.out.println("Our algorithm took " + dca + " milliseconds");

		if (res.feasible) {
			System.out.println("yes, the problem is feasible");
			//System.out.println(Arrays.toString(res.sol));
		} else {
			System.out.println("no, the problem is infeasible");
		}
		*/

		

		/********************************************************************************
        **
        **MDA: solve the problem by FastMDA and record the time
        **
        ********************************************************************************/
      	//
		RAPNC testFastMDA = new RAPNC(obj, lbVar, ubVar, lbNested, ubNested);
		startTime = System.currentTimeMillis();
		ResultTypeMDA resFastMDA = testFastMDA.FastMDA();
		//System.out.println("MDA Done!");
		//ResultTypeMDA resFastMDA = new ResultTypeMDA(resFlow.sol, resFlow.sol, resFlow.sol, resFlow.sol);

		//time 
		endTime = System.currentTimeMillis();
		long timeFastMDA = endTime - startTime;

		//number of subproblems
		long number_subproblem_MDA = testFastMDA.number_subproblem;


    	//check whether the two algorithm produce the same solution
    	for (int i = 0; i < dimension; i++) {
			if (Math.abs(res.sol[i] - resFastMDA.aa[i]) >= 0.01) {
				System.out.println("No, the solution x[" + i + "] is different.");
				System.out.println(res.sol[i] + " " + " " + resFastMDA.aa[i]);
				System.out.println("Coefficient: " + obj.get(i).getValue(3) + obj.get(i).getValue(2));
			}
			
		}
	
		System.out.println(
				String.format("%10s", dimension) 
				+ String.format("%10s", ((double) timeDCA) / 1000) 
				+ String.format("%10s", ((double) timeFastMDA) / 1000) 
				+ String.format("%20s", (number_subproblem_DCA)) 
				+ String.format("%20s", (number_subproblem_MDA)) 
		);
    }
}
