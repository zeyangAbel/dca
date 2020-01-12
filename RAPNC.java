/*******************************************************************************
 *  Compilation:  javac -cp ./ RAPNC.java
 *  Execution:    java -cp ./ RAPNC
 *  
 *  A class for resource allocation problem with lower and upper nested constraints (DRAP-NC) 
 *  @author Zeyang Wu @ University of Minnesota
 *
 *******************************************************************************/

/******************************************************************************
* This is a class for RAPNC with separable convex objectives
* It provides three methods to solve the problem      
* Method: 
*     1. solveInteger()-----DCA: (divide and conquer algorithm) 
*     This algorithm is inspired by He et al. 2017 paper Speed optimization over a path with heterogeneous arc costs.
*     The algorithm solve a RAP relaxation of RAP-NC and selects the maximum violation of the nested constraints
*     Then the problem can be divide into two subproblems in the form of RAP-NC and are solved recursively.
*     Running time: Theta(n^2 log B)
*     Benefit: easy implementation and no Lipschitz continuity is required
*     solveIntegerLinear() ---- With the O(n) subroutine for solving DRAP, the algorithm can be sped up to Theta(n^2)
*     @param RAP class must be compiled
*
*     2. FastMDA(int u, int v)------MDA: (monotonic decomposition algorithm)
*     This algorithm is proposed by Vidal et al's 2018 paper on IJOC: 
*     Separable Convex Optimization with Nested Lower and Upper Constraints
*     The original implementation is MDA(int u, int v) The main issue is that it requires an Lipschitz constant
*     After a discussion on INFROMS 2018, we learn that the Lipschitz constant is not mandatory and the algorithm can be improved by discussing 
*     the boundary cases. Hence, we improve it to and implement FastMDA() 
*     Running time: O(n log n log B)
*     LinearMDA() ---- With the O(n) subroutine for solving DRAP, the algorithm can be sped up to O(n log n)
*     
*     @param RAP class must be compiled
*         
*/

import java.util.*;

public class RAPNC {
	//problem dimension
	int dimension;

    //list of function oracles 
	List<Function> obj;

	//list of parameters
	long[] ubVar;
	long[] lbVar;
	long[] ubNested;
	long[] lbNested;
    long[] lbCopyMDA; //reserve bounds for MDA
    long[] ubCopyMDA;
	private long scaleFactor;
    long number_subproblem;

    //normal constructor
    public RAPNC(int K) {
    	//normal constructor
    	this.obj = new ArrayList<Function>();
		this.lbVar = new long[K];
		this.ubVar =  new long[K];
		this.lbNested =  new long[K];
		this.ubNested =  new long[K];
		this.dimension = lbVar.length;
		this.scaleFactor = 1;
        this.number_subproblem = 0;
    }

    //normal constructor
	public RAPNC(List<Function> obj, long[] lbVar, long[] ubVar, long[] lbNested, long[] ubNested) {
		this.obj = obj;
		this.lbVar = lbVar;
		this.ubVar = ubVar;
		this.lbNested = lbNested;
		this.ubNested = ubNested;
		this.dimension = lbVar.length;
		this.scaleFactor = 1;		
        this.number_subproblem = 0;
	}

    //This is used in solving continuous RAP-NC    
    public void setScaleFactor(long scaleFactor) {
    	this.scaleFactor = scaleFactor;
    }

     /**
     * createRAP()
     * By this method we can create an instance of RAP by relaxing the nested constraints.
     * Time-Complexity: O(n)     *
     * @param no argument
     * @return RAP relaxation of RAP-NC
     */
    public RAP createRAP() {
    	//The resource bound is the bound of the last nested constraint.
    	long rapB = ubNested[dimension - 1];
    	long[] raplb = lbVar.clone();
    	long[] rapub = ubVar.clone();
    	RAP res = new RAP(obj, rapB, raplb, rapub);
    	res.scaleFactor = this.scaleFactor;
    	return res;
    }

    /**
     * createRAPNC
     * By this method we can create two RAPNC subproblems by a index K.
     * RAP-NC(s, K) and RAP-NC(K + 1, e) 
     * Time-Complexity: O(n) 
     * 
     * @return List<RAPNC> containing the two subproblems
     */
    public List<RAPNC> createRAPNC(int K) {
    	//The first subproblem contains K + 1 variables and the second subproblem contains (n - K - 1) variables
    	//return a list containing the two problems

    	//copy array
    	RAPNC left = new RAPNC(K + 1);
    	RAPNC right = new RAPNC(dimension - K - 1);
    	//setup left obj
    	System.arraycopy(this.ubNested, 0, left.ubNested, 0, K + 1);
    	System.arraycopy(this.lbNested, 0, left.lbNested, 0, K + 1);
    	System.arraycopy(this.ubVar, 0, left.ubVar, 0, K + 1);
    	System.arraycopy(this.lbVar, 0, left.lbVar, 0, K + 1);
    	left.dimension = K + 1;
    	left.obj = this.obj;
    	left.scaleFactor = this.scaleFactor;


    	//setup right
    	right.dimension = dimension - K - 1;
    	System.arraycopy(this.ubVar, K + 1, right.ubVar, 0, dimension - K - 1);
    	System.arraycopy(this.lbVar, K + 1, right.lbVar, 0, dimension - K - 1);
    	right.obj = new ArrayList<Function>();
    	for (int i = K + 1; i < dimension; i++) {
    		right.ubNested[i - K - 1] = ubNested[i] - ubNested[K];
    		right.lbNested[i - K - 1] = lbNested[i] - ubNested[K];
    		right.obj.add(this.obj.get(i));
    	}
    	right.scaleFactor = this.scaleFactor;

    	List<RAPNC> res = new ArrayList<RAPNC>();
    	res.add(left);
    	res.add(right);
    	return res;
    }


    /**
     * solveIntegerDCA() 
     * This is the implementation of DCA method
     * Time-Complexity: O(n^2 log(B)) 
     * 
     * @return ResultTypeRAPNC containing the solution and feasibility
     */
    public ResultTypeRAPNC solveIntegerDCA() {
    	//Solve the RAPNC with integer variables
    	//1. solve the relaxation problem
    	//2. find the maximum violation and then divide the problem into two subproblems
    	//3. solve the two subproblem recursively and conquer the results

        this.number_subproblem++; //record the subproblems
    	//Trivial case
    	if (dimension == 1) {
    		//check feasibility
    		if (ubNested[0] >= lbVar[0] && ubNested[0] <= ubVar[0]) {
    			long[] sol = new long[1];
    			sol[0] = ubNested[0];
    			return new ResultTypeRAPNC(true, sol);
    		} else {
    			return new ResultTypeRAPNC(false, null);
    		}
    	}


    	//call Greedy Algorithm to solve the relaxed instance with greedy algorithm.
    	RAP re = createRAP();
    	ResultTypeRAP solRAP = re.solveRAP(); 

    	if (!solRAP.feasible) {
    		return new ResultTypeRAPNC(false, null);
    	}
    	//obtain a solution to RAP
    	long[] solRe = solRAP.sol.clone();
    	/*Debug
        
		System.out.println("yes, the problem is feasible");
		System.out.println(Arrays.toString(solRe));
		*/

    	//find the maximum violation
    	long sum = 0;
    	int maxIndex = -1;
    	long maxVio = 0;
    	int maxFlag = 0;//excess 1, shortage 0
    	long[] violation = new long[dimension];
    	
    	for (int i = 0; i < dimension; i++) {
    		sum += solRe[i];
    		int flag = 1;
    		if (sum > ubNested[i]) {
    			violation[i] = sum - ubNested[i];
    		} else if (sum < lbNested[i]) {
    			violation[i] = lbNested[i] - sum;
    			flag = 0;
    		}
    		if (violation[i] > maxVio) {
    			maxIndex = i;
    			maxVio = violation[i];
    			maxFlag = flag;
    		}
    	}

        //If the solution to RAP satisfies all nested constraints. 
    	if (maxIndex == -1) {
    		//.out.println(Arrays.toString(solRe));
    		return new ResultTypeRAPNC(true, solRe);
    	} 
    	/*Debug
    	System.out.print("Split at:");
        System.out.print(maxIndex);
        System.out.print(" with ");
        System.out.println(maxFlag);
		*/

    	//else divide the problem into two small problems. 
    	if (maxFlag == 0) {
    		ubNested[maxIndex] = lbNested[maxIndex];
    	} else {
    		lbNested[maxIndex] = ubNested[maxIndex];
    	}
    	
    	/*Debug
    	System.out.print("TightenNested:");
        System.out.print(ubNested[maxIndex]);
        System.out.print(" with ");
        System.out.println(ubNested[maxIndex]);
    	*/  	
    	
    	List<RAPNC> divide = createRAPNC(maxIndex);
    	ResultTypeRAPNC left = divide.get(0).solveIntegerDCA(); 
    	ResultTypeRAPNC right = divide.get(1).solveIntegerDCA();

        this.number_subproblem += divide.get(0).number_subproblem;
        this.number_subproblem += divide.get(1).number_subproblem;

    	if (!(left.feasible && right.feasible)) {
    		return new ResultTypeRAPNC(false, null);
    	} else {
    		System.arraycopy(left.sol, 0, solRe, 0, left.sol.length);
    		System.arraycopy(right.sol, 0, solRe, maxIndex + 1, right.sol.length);
    		return new ResultTypeRAPNC(true, solRe);
    	}
    } 

    /**
     * solveIntegerLinearDCA() 
     * This is the implementation of DCA method for linear objectives
     * Time-Complexity: O(n^2) 
     * 
     * @return ResultTypeRAPNC containing the solution and feasibility
     */
    public ResultTypeRAPNC solveIntegerLinearDCA() {
        //Solve the RAPNC with integer variables
        //1. solve the relaxation problem
        //2. find the maximum violation and then divide the problem into two subproblems
        //3. solve the two subproblem recursively and conquer the results

        //Trivial case
        if (dimension == 1) {
            //check feasibility
            if (ubNested[0] >= lbVar[0] && ubNested[0] <= ubVar[0]) {
                long[] sol = new long[1];
                sol[0] = ubNested[0];
                return new ResultTypeRAPNC(true, sol);
            } else {
                return new ResultTypeRAPNC(false, null);
            }
        }


        //call Greedy Algorithm to solve the relaxed instance with greedy algorithm.
        RAP re = createRAP();
        ResultTypeRAP solRAP = re.solveRAPLinear(); 

        if (!solRAP.feasible) {
            return new ResultTypeRAPNC(false, null);
        }
        //obtain a solution to RAP
        long[] solRe = solRAP.sol.clone();
        /*Debug
        System.out.println("yes, the problem is feasible");
        System.out.println(Arrays.toString(solRe));
        */

        //find the maximum violation
        long sum = 0;
        int maxIndex = -1;
        long maxVio = 0;
        int maxFlag = 0;//excess 1, shortage 0
        long[] violation = new long[dimension];
        
        for (int i = 0; i < dimension; i++) {
            sum += solRe[i];
            int flag = 1;
            if (sum > ubNested[i]) {
                violation[i] = sum - ubNested[i];
            } else if (sum < lbNested[i]) {
                violation[i] = lbNested[i] - sum;
                flag = 0;
            }
            if (violation[i] > maxVio) {
                maxIndex = i;
                maxVio = violation[i];
                maxFlag = flag;
            }
        }

        //If the solution to RAP satisfies all nested constraints. 
        if (maxIndex == -1) {
            //System.out.println(Arrays.toString(solRe));
            return new ResultTypeRAPNC(true, solRe);
        } 
        /*Debug
        System.out.print("Split at:");
        System.out.print(maxIndex);
        System.out.print(" with ");
        System.out.println(maxFlag);
        */

        //else divide the problem into two small problems. 
        if (maxFlag == 0) {
            ubNested[maxIndex] = lbNested[maxIndex];
        } else {
            lbNested[maxIndex] = ubNested[maxIndex];
        }
        
        /*Debug
        System.out.print("TightenNested:");
        System.out.print(ubNested[maxIndex]);
        System.out.print(" with ");
        System.out.println(ubNested[maxIndex]);
        */      
        
        List<RAPNC> divide = createRAPNC(maxIndex);
        ResultTypeRAPNC left = divide.get(0).solveIntegerLinearDCA(); 
        ResultTypeRAPNC right = divide.get(1).solveIntegerLinearDCA();

        if (!(left.feasible && right.feasible)) {
            return new ResultTypeRAPNC(false, null);
        } else {
            System.arraycopy(left.sol, 0, solRe, 0, left.sol.length);
            System.arraycopy(right.sol, 0, solRe, maxIndex + 1, right.sol.length);
            return new ResultTypeRAPNC(true, solRe);
        }
    }  

   


    /**
     * MDA() 
     * This is the implementation of MDA method，it initiate the Main recursion MDA(u,v)
     * Time-Complexity: O(n log(n) log(B)) 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @param obj the function oracles should be transformed to corresponding obj with penalties
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */
    public ResultTypeMDA MDA() {
        lbCopyMDA = Arrays.copyOfRange(lbVar, 0, dimension);
        ubCopyMDA = Arrays.copyOfRange(ubVar, 0, dimension);
        return MDA(0, dimension - 1);
    }

    /**
     * MDA() 
     * This is the implementation of MDA method
     * Time-Complexity: O(n log(n) log(B)) 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @param obj the function oracles should be transformed to corresponding obj with penalties
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @param v start index
     * @param u end index
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */

    public ResultTypeMDA MDA(int v, int w) {
        
        //initialize return of the f
        long[] aa = new long[w - v + 1];
        long[] ab = new long[w - v + 1];
        long[] ba = new long[w - v + 1];
        long[] bb = new long[w - v + 1];

        //simple cases
        if (v == w) {
            if (v == 0) {
                aa[0] = lbNested[v];
                ab[0] = ubNested[v];
                ba[0] = lbNested[v];
                bb[0] = ubNested[v];
            } else {
                aa[0] = lbNested[v] - lbNested[v - 1];
                ab[0] = ubNested[v] - lbNested[v - 1];
                ba[0] = lbNested[v] - ubNested[v - 1];
                bb[0] = ubNested[v] - ubNested[v - 1];
            }
            return new ResultTypeMDA(aa, ab, ba, bb);            
        }

        //divide 
        int u = v + (w - v) / 2;
        ResultTypeMDA left = MDA(v, u);
        ResultTypeMDA right = MDA(u + 1, w);

        //conquer
        //update the bounds
        //aa
        //check the conditions
        boolean flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        boolean flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            aa = subproblemRAPSolve(v, w, lbNested[w]);
        } else {
            aa = subproblemRAPSolve(v, w, lbNested[w] - lbNested[v - 1]);
        }

        //avoid the last computation
        if (v == 0 && w == dimension - 1) {
            //System.out.println("yes");
            return new ResultTypeMDA(aa, aa, aa, aa);
        }

        //ab
        flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            ab = subproblemRAPSolve(v, w, ubNested[w]);
        } else {
            ab = subproblemRAPSolve(v, w, ubNested[w] - lbNested[v - 1]);
        }

        //ba
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            ba = subproblemRAPSolve(v, w, lbNested[w]);
        } else {
            ba = subproblemRAPSolve(v, w, lbNested[w] - ubNested[v - 1]);
        }

        //bb
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            bb = subproblemRAPSolve(v, w, ubNested[w]);
        } else {
            bb = subproblemRAPSolve(v, w, ubNested[w] - ubNested[v - 1]);
        }


        //return 
        return new ResultTypeMDA(aa, ab, ba, bb); 
    }

    /**
     * LinearMDA() 
     * This is for RAP-NC with linear objectives
     * This is the implementation of MDA method，it initiate the Main recursion LinearMDA(u,v)
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * @param obj the function oracles should be transformed to corresponding obj with penalties
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * Time-Complexity: O(n log(n)) 
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */
    public ResultTypeMDA LinearMDA() {
        lbCopyMDA = Arrays.copyOfRange(lbVar, 0, dimension);
        ubCopyMDA = Arrays.copyOfRange(ubVar, 0, dimension);
        return LinearMDA(0, dimension - 1);
    }

    /**
     * LinearMDA() 
     * This is the implementation of LinearMDA method. 
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1
     * The @code of this method is incorrect.
     * !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     * Time-Complexity: O(n log(n)) 
     * @param v start index
     * @param u end index
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */

    public ResultTypeMDA LinearMDA(int v, int w) {
        
        //initialize return of the solution
        long[] aa = new long[w - v + 1];
        long[] ab = new long[w - v + 1];
        long[] ba = new long[w - v + 1];
        long[] bb = new long[w - v + 1];

        //simple cases
        if (v == w) {
            if (v == 0) {
                aa[0] = lbNested[v];
                ab[0] = ubNested[v];
                ba[0] = lbNested[v];
                bb[0] = ubNested[v];
            } else {
                aa[0] = lbNested[v] - lbNested[v - 1];
                ab[0] = ubNested[v] - lbNested[v - 1];
                ba[0] = lbNested[v] - ubNested[v - 1];
                bb[0] = ubNested[v] - ubNested[v - 1];
            }
            return new ResultTypeMDA(aa, ab, ba, bb);            
        }

        //divide 
        int u = v + (w - v) / 2;
        ResultTypeMDA left = LinearMDA(v, u);
        ResultTypeMDA right = LinearMDA(u + 1, w);

        //conquer
        //update the bounds
        //aa
        //check the conditions
        boolean flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        boolean flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            aa = subproblemRAPSolveLinear(v, w, lbNested[w]);
        } else {
            aa = subproblemRAPSolveLinear(v, w, lbNested[w] - lbNested[v - 1]);
        }

        //avoid the last computation
        if (v == 0 && w == dimension - 1) {
            //System.out.println("yes");
            return new ResultTypeMDA(aa, aa, aa, aa);
        }

        //ab
        flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            ab = subproblemRAPSolveLinear(v, w, ubNested[w]);
        } else {
            ab = subproblemRAPSolveLinear(v, w, ubNested[w] - lbNested[v - 1]);
        }

        //ba
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            ba = subproblemRAPSolveLinear(v, w, lbNested[w]);
        } else {
            ba = subproblemRAPSolveLinear(v, w, lbNested[w] - ubNested[v - 1]);
        }

        //bb
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            bb = subproblemRAPSolveLinear(v, w, ubNested[w]);
            long[] bstar = subproblemRAPSolve(v, w, ubNested[w]);

            for (int i = 0; i < bb.length; i++) {
                if ((Math.abs(bb[i] - bstar[i]) >= 0.01)) {
                System.out.println("No, the solution x[" + i + "] is different.");
                System.out.print(obj.get(i).getValue(1));
                System.out.print(" ");
                System.out.print(lbCopyMDA[i]);
                System.out.print(" ");
                System.out.print(ubCopyMDA[i]);
                System.out.print(" ");
                System.out.print(bb[i]);
                System.out.print(" ");
                System.out.println(bstar[i]);     
                }
            }
        } else {
            bb = subproblemRAPSolveLinear(v, w, ubNested[w] - ubNested[v - 1]);
        }


        //return 
        return new ResultTypeMDA(aa, ab, ba, bb); 
    }

    /**
     * check method 
     * This is an method check if the component-wise vector inequality nums1 <= nums2 is violated 
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return boolean variable if nums1 <= nums2 then FALSE
     */
    public boolean check(long[] nums1, long[] nums2, int v, int u) {
        //check if nums1 <= nums2
        boolean res = false;
        for (int i = v; i <= u; i++) {
            if (nums1[i - v] > nums2[i - v]) {
                res = true;
            }
        }
        return res;
    }

    /**
     * adjust method 
     * This is an method adjust the two vector such that nums1 <= nums2 is satisfied in a way that
     * the modified arrays are still optimal to the corresponding subproblems
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return no return 
     */
    public void adjust(long[] nums1, long[] nums2, int v, int u) {
        long delta = 0;
        for (int i = v; i <= u; i++) {
             if (nums1[i - v] > nums2[i - v]) {
                delta += nums1[i - v] - nums2[i - v];
                nums1[i - v] = nums2[i - v];
             }
        }

        for (int i = v; i <= u; i++) {
             if (nums1[i - v] < nums2[i - v]) {
                long ins = Math.min(nums2[i - v] - nums1[i - v], delta);
                delta -= ins;
                nums1[i - v] += ins;
             }
        }
    }

    /**
     * update method 
     * This is an method used to update the variable bounds in the MDA method
     * Time-Complexity: O(n) 
     * @param nums1 
     * @param nums2 
     * @param v start index
     * @param u end index
     * @return no
     */
    public void updateBounds(long[] nums1, long[] nums2, int v, int u) {
        for (int i = v; i <= u; i++) {
             lbCopyMDA[i] = nums1[i - v];
             ubCopyMDA[i] = nums2[i - v];
        }
    }

    /**
     * subproblemRAPSolve 
     * This is an method used to solve the special RAP subproblems in MDA method
     * Time-Complexity: O(n) 
     * @param v start index
     * @param u end index
     * @param LR the total amount of resource in the subproblems
     * @return solution to the subproblem. Here we simply assume that the problem is always feasible
     */
    public long[] subproblemRAPSolve(int v, int w, long LR) { 
        long rapB = LR;
        long[] raplb = Arrays.copyOfRange(lbCopyMDA, v, w + 1);
        long[] rapub = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
        List<Function> rapObj = obj.subList(v, w + 1);
        RAP rap = new RAP(rapObj, rapB, raplb, rapub);
        rap.scaleFactor = this.scaleFactor;
        ResultTypeRAP res = rap.solveRAP();
        if (!res.feasible) {
            System.out.println("Subproblem" + v + " " + w  + "Infeasible");
        }
        return res.sol;
    }

    /**
     * subproblemRAPSolveLinear
     * This is an method used to solve the special RAP subproblems in MDA method
     * Time-Complexity: O(n) 
     * @param v start index
     * @param u end index
     * @param LR the total amount of resource in the subproblems
     * @return solution to the subproblem. Here we simply assume that the problem is always feasible
     */
    public long[] subproblemRAPSolveLinear(int v, int w, long LR) { 
        long rapB = LR;

        long sumUb = 0;
        long sumLb = 0;
        long diffUb = 0;
        long diffLb = 0;

        long[] raplb = Arrays.copyOfRange(lbCopyMDA, v, w + 1);
        long[] rapub = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
        
        //This is the implementation of the formula in Vidal's paper
        for (int i = v; i < w + 1; i++) {
            sumUb += ubVar[i];
            sumLb += lbVar[i];
            diffUb += rapub[i - v] - ubVar[i];
            diffLb += raplb[i - v] - lbVar[i];
            if (raplb[i - v] < lbVar[i]) {
                raplb[i - v] = lbVar[i];
            }
            if (rapub[i - v] > ubVar[i]) {
                rapub[i - v] = ubVar[i];
            }
        }

        ResultTypeRAP res = new ResultTypeRAP(true, new long[w - v + 1]);
        if (sumUb < LR) {
            long ratio = (LR - sumUb) / (diffUb);
            for (int i = v; i < w + 1; i++) {
                res.sol[i - v] = ubVar[i] - ratio * (ubCopyMDA[i] - ubVar[i]);
            }
        } else if (sumLb > LR) {
            long ratio = (LR - sumLb) / (diffLb);
            for (int i = v; i < w + 1; i++) {
                res.sol[i - v] = lbVar[i] - ratio * (lbCopyMDA[i] - lbVar[i]);
            }
        } else {
            List<Function> rapObj = obj.subList(v, w + 1);
            RAP rap = new RAP(rapObj, rapB, raplb, rapub);
            rap.scaleFactor = this.scaleFactor;
            res = rap.solveRAPLinear();
            if (!res.feasible) {
                System.out.println("Subproblem" + v + " " + w  + "Infeasible");
            }
        }
       
        return res.sol;
    }

    // New implementation of MDA 
    /**
     * FastMDA() 
     * This is an implementation of MDA method，it initiate the Main recursion FastMDA(u,v)
     * We discuss with vidal at informs annual meeting 2018 to have a better understanding of the detailed implementation.
     * Time-Complexity: O(n log(n) log(B)) 
     * @param obj the function oracles should be transformed to corresponding obj with penalties
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */
    public ResultTypeMDA FastMDA() {
        //lbCopyMDA = Arrays.copyOfRange(lbVar, 0, dimension);
        //ubCopyMDA = Arrays.copyOfRange(ubVar, 0, dimension);
        lbCopyMDA = new long[dimension];
        ubCopyMDA = new long[dimension];
        //Arrays.fill(lbCopyMDA, - 2 * ubNested[dimension - 1]);
        //Arrays.fill(ubCopyMDA, 2 * ubNested[dimension - 1]);
        return FastMDA(0, dimension - 1);
    }

    /**
     * FastMDA() 
     * This is the implementation of MDA method
     * Time-Complexity: O(n log(n) log(B)) 
     * @param obj 
     * @param v start index
     * @param u end index
     * @return ResultTypeMDA containing the solution to the four subproblems and feasibility
     */

    public ResultTypeMDA FastMDA(int v, int w) {
        
        this.number_subproblem += 4; // record the number of subproblems
        //initialize return of the f
        long[] aa = new long[w - v + 1];
        long[] ab = new long[w - v + 1];
        long[] ba = new long[w - v + 1];
        long[] bb = new long[w - v + 1];

        //simple cases
        if (v == w) {
            if (v == 0) {
                aa[0] = lbNested[v];
                ab[0] = ubNested[v];
                ba[0] = lbNested[v];
                bb[0] = ubNested[v];
            } else {
                aa[0] = lbNested[v] - lbNested[v - 1];
                ab[0] = ubNested[v] - lbNested[v - 1];
                ba[0] = lbNested[v] - ubNested[v - 1];
                bb[0] = ubNested[v] - ubNested[v - 1];
            }
            return new ResultTypeMDA(aa, ab, ba, bb);            
        }

        /*
        An new idea to improve the algorithm: We solve the subproblems before applying a divide step. However it is slower.
        */
        /*
        //solve the subproblems before applying a divide step
        if (v == 0) {
            aa = subproblemRAPSolveFastMDAAddChecker(v, w, lbNested[w]);
            ab = subproblemRAPSolveFastMDAAddChecker(v, w, ubNested[w]);
            ba = subproblemRAPSolveFastMDAAddChecker(v, w, lbNested[w]);
            bb = subproblemRAPSolveFastMDAAddChecker(v, w, ubNested[w]);
        } else {
            aa = subproblemRAPSolveFastMDAAddChecker(v, w, lbNested[w] - lbNested[v - 1]);
            ab = subproblemRAPSolveFastMDAAddChecker(v, w, ubNested[w] - lbNested[v - 1]);
            ba = subproblemRAPSolveFastMDAAddChecker(v, w, lbNested[w] - ubNested[v - 1]);
            bb = subproblemRAPSolveFastMDAAddChecker(v, w, ubNested[w] - ubNested[v - 1]);
            //System.out.println("Something may be wrong");
        }

        //check if there is any violations
        boolean has_violation = false;  

        long aa_sum = 0;
        long ab_sum = 0;
        long ba_sum = 0;
        long bb_sum = 0;

        long la_adjust = 0;
        long lb_adjust = 0;
        if (v != 0) {
            la_adjust = lbNested[v - 1];
            lb_adjust = ubNested[v - 1];
        }

        //under original bounds, the problem may not be feasible
        if (aa == null || ab == null || ba == null || bb == null) {
            has_violation = true;
        }

        for (int i = v; i <= w; i++) {

            if (has_violation) {
                break;
            }

            aa_sum += aa[i - v];
            ab_sum += ab[i - v];
            ba_sum += ba[i - v];
            bb_sum += bb[i - v];

            if (!(aa_sum >= lbNested[i] - la_adjust && aa_sum <= ubNested[i] - la_adjust)) {
                has_violation = true;
                break;
            }   
            if (!(ab_sum >= lbNested[i] - la_adjust && ab_sum <= ubNested[i] - la_adjust)) {
                has_violation = true;
                break;
            }
            if (!(ba_sum >= lbNested[i] - lb_adjust && ba_sum <= ubNested[i] - lb_adjust)) {
                has_violation = true;
                break;
            }
            if (!(bb_sum >= lbNested[i] - lb_adjust && bb_sum <= ubNested[i] - lb_adjust)) {
                has_violation = true;
                break;
            }
            
        }


        if (!has_violation) {
            return new ResultTypeMDA(aa, ab, ba, bb);
        }

        this.number_subproblem += 4; // record the number of subproblems
        */
        /*
        New idea to test.
        */



        //divide 
        int u = v + (w - v) / 2;
        ResultTypeMDA left = FastMDA(v, u);
        ResultTypeMDA right = FastMDA(u + 1, w);


        //conquer
        //update the bounds
        //aa
        //check the conditions
        boolean flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        boolean flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            aa = subproblemRAPSolveFastMDA(v, w, lbNested[w]);
        } else {
            aa = subproblemRAPSolveFastMDA(v, w, lbNested[w] - lbNested[v - 1]);
        }

        //avoid the last computation
        if (v == 0 && w == dimension - 1) {
            //System.out.println("yes");
            return new ResultTypeMDA(aa, aa, aa, aa);
        }

        //ab
        flag = check(left.aa, left.ab, v, u);
        if (flag) {
            adjust(left.aa, left.ab, v, u);
        }
        updateBounds(left.aa, left.ab, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            ab = subproblemRAPSolveFastMDA(v, w, ubNested[w]);
        } else {
            ab = subproblemRAPSolveFastMDA(v, w, ubNested[w] - lbNested[v - 1]);
        }

        //ba
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.ba, right.aa, u + 1, w);
        if (flag2) {
            adjust(right.ba, right.aa, u + 1, w);
        }
        updateBounds(right.ba, right.aa, u + 1, w);
        if (v == 0) {
            ba = subproblemRAPSolveFastMDA(v, w, lbNested[w]);
        } else {
            ba = subproblemRAPSolveFastMDA(v, w, lbNested[w] - ubNested[v - 1]);
        }

        //bb
        //check the conditions
        flag = check(left.ba, left.bb, v, u);
        if (flag) {
            adjust(left.ba, left.bb, v, u);
        }
        updateBounds(left.ba, left.bb, v, u);

        flag2 = check(right.bb, right.ab, u + 1, w);
        if (flag2) {
            adjust(right.bb, right.ab, u + 1, w);
        }
        updateBounds(right.bb, right.ab, u + 1, w);
        if (v == 0) {
            bb = subproblemRAPSolveFastMDA(v, w, ubNested[w]);
        } else {
            bb = subproblemRAPSolveFastMDA(v, w, ubNested[w] - ubNested[v - 1]);
        }


        //return 
        return new ResultTypeMDA(aa, ab, ba, bb); 
    }

    /**
     * subproblemRAPSolveFastMDA2
     * This is an method used to solve the special RAP subproblems in MDA method
     * Time-Complexity: O(n log n) 
     * @param v start index
     * @param u end index
     * @param LR the total amount of resource in the subproblems
     * @return solution to the subproblem. Here we simply assume that the problem is always feasible
     */
    public long[] subproblemRAPSolveFastMDA2(int v, int w, long LR) { 
        //set up checkers
        long sum_c_bar = 0;
        long sum_c = 0;
        long sum_d_bar = 0;
        long sum_d = 0;

        long rapB = LR;
        long[] raplb = new long[w + 1 - v];
        long[] rapub = new long[w + 1 - v];
        List<Function> rapObj = obj.subList(v, w + 1);

        for (int i = v; i < w + 1; i++) {
            sum_c_bar += lbCopyMDA[i];
            sum_c += lbVar[i];
            sum_d_bar += ubCopyMDA[i];
            sum_d += ubVar[i];

            if (lbCopyMDA[i] > ubCopyMDA[i]) {
                System.out.println("Something of c_bar and d_bar may be wrong");
                //here is correct.
            }
        }

        if (sum_c_bar > LR || sum_d_bar < LR) {
            System.out.println("Weired"); // These constraints are always fulfilled.
        }

        //case one: the problem is infeasible w.r.t. sum_c
        // This can be viewed as a special implementation of the greedy algorithms for RAP. (Because the penalty cost are all the same and linear.)
        if (sum_c_bar <= LR && sum_c > LR) {
            // return some solution
            long[] construct_sol = Arrays.copyOfRange(lbCopyMDA, v, w + 1);

            for (int i = v; i < w + 1; i++) {
                
                if (LR - sum_c_bar == 0) {
                    break;
                }

                long delta = Math.max(Math.min(lbVar[i] - lbCopyMDA[i], LR - sum_c_bar), 0);
                delta = Math.min(ubCopyMDA[i] - lbCopyMDA[i], delta);
                construct_sol[i - v] += delta;
                sum_c_bar += delta;
            }
            if (LR - sum_c_bar > 0) {
                System.out.println("Not a feasible solution."); 
            }
            return construct_sol;
        }

        //case two:
        if (sum_d_bar >= LR && sum_d < LR) {
            // return some solution
            long[] construct_sol = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
            //System.out.println("Case 2 Checked");
            for (int i = v; i < w + 1; i++) {
                
                if (sum_d_bar - LR == 0) {
                    break;
                }
                
                long delta = Math.max(Math.min(ubCopyMDA[i] - ubVar[i], sum_d_bar - LR), 0);
                delta = Math.min(ubCopyMDA[i] - lbCopyMDA[i], delta);
                construct_sol[i - v] -= delta;
                sum_d_bar -= delta;
            }

            if (sum_d_bar - LR > 0) {
                System.out.println("Not a feasible solution.");
            }
            return construct_sol;
        }


        //case three: the RAP problem is feasible under the original bound
        //raplb = Arrays.copyOfRange(lbVar, v, w + 1);
        //rapub = Arrays.copyOfRange(ubVar, v, w + 1);
        raplb = Arrays.copyOfRange(lbCopyMDA, v, w + 1);
        rapub = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
        RAP rap = new RAP(rapObj, rapB, raplb, rapub);
        sum_c = 0;
        sum_d = 0;
        for (int i = v; i < w + 1; i++) {
            sum_c += raplb[i - v];
            sum_d += rapub[i - v];
            if (raplb[i - v] > rapub[i - v]) {
                System.out.println("Something may be wrong");
            }
        }
        if (sum_c > rapB) {
            System.out.println("lower bound is wrong");
        } else if (sum_d < rapB) {
            System.out.println("upper bound is wrong");
        }
        //System.out.println("Checked");
        rap.scaleFactor = this.scaleFactor;
        ResultTypeRAP res = rap.solveRAP();
        if (!res.feasible) {
            System.out.println("Subproblem" + v + " " + w + "Infeasible");
        }
        return res.sol;
    }


    /**
     * subproblemRAPSolveFastMDA
     * This is an method used to solve the special RAP subproblems in MDA method
     * Time-Complexity: O(n log n) 
     * @param v start index
     * @param u end index
     * @param LR the total amount of resource in the subproblems
     * @return solution to the subproblem. Here we simply assume that the problem is always feasible
     */
    public long[] subproblemRAPSolveFastMDA(int v, int w, long LR) { 
        //set up checkers
        long sum_c_prime = 0;
        long sum_d_prime = 0;
        long sum_c_bar = 0;
        long sum_d_bar = 0;

        long rapB = LR;
        long[] raplb = new long[w + 1 - v];
        long[] rapub = new long[w + 1 - v];
        List<Function> rapObj = obj.subList(v, w + 1);

        for (int i = v; i < w + 1; i++) {

            raplb[i - v] = Math.max(lbCopyMDA[i], Math.min(lbVar[i], ubCopyMDA[i]));
            rapub[i - v] = Math.min(ubCopyMDA[i], Math.max(ubVar[i], lbCopyMDA[i]));
            
            sum_c_bar += lbCopyMDA[i];
            sum_d_bar += ubCopyMDA[i];
            sum_c_prime += raplb[i - v];      
            sum_d_prime += rapub[i - v];
        }

        //case one: the problem is infeasible w.r.t. sum_c
        // This can be viewed as a special implementation of the greedy algorithms for RAP. (Because the penalty cost are all the same and linear.)
        if (sum_c_prime > LR) {
            // return some solution
            long[] construct_sol = Arrays.copyOfRange(lbCopyMDA, v, w + 1);

            for (int i = v; i < w + 1; i++) {
                
                if (LR - sum_c_bar == 0) {
                    break;
                }

                long delta = Math.max(Math.min(raplb[i - v] - lbCopyMDA[i], LR - sum_c_bar), 0);
                construct_sol[i - v] += delta;
                sum_c_bar += delta;
            }

            return construct_sol;
        }

        //case two:
        if (sum_d_prime < LR) {
            // return some solution
            long[] construct_sol = Arrays.copyOfRange(ubCopyMDA, v, w + 1);
            //System.out.println("Case 2 Checked");
            for (int i = v; i < w + 1; i++) {
                
                if (sum_d_bar - LR == 0) {
                    break;
                }
                
                long delta = Math.max(Math.min(ubCopyMDA[i] - rapub[i - v], sum_d_bar - LR), 0);
                construct_sol[i - v] -= delta;
                sum_d_bar -= delta;
            }

            return construct_sol;
        }


        //case three: the RAP problem is feasible under the original bound
        RAP rap = new RAP(rapObj, rapB, raplb, rapub);
        long sum_c = 0;
        long sum_d = 0;
        for (int i = v; i < w + 1; i++) {
            sum_c += raplb[i - v];
            sum_d += rapub[i - v];
            if (raplb[i - v] > rapub[i - v]) {
                System.out.println("Something may be wrong");
            }
        }
        if (sum_c > rapB) {
            System.out.println("lower bound is wrong");
        } else if (sum_d < rapB) {
            System.out.println("upper bound is wrong");
        }
        //System.out.println("Checked");
        rap.scaleFactor = this.scaleFactor;
        ResultTypeRAP res = rap.solveRAP();
        //ResultTypeRAP res = rap.solveRAPLinear();
        if (!res.feasible) {
            System.out.println("Subproblem" + v + " " + w + "Infeasible");
        }
        return res.sol;
    }

    public long[] subproblemRAPSolveFastMDAAddChecker(int v, int w, long LR) { 
        //set up checkers
        long sum_c_prime = 0;
        long sum_d_prime = 0;
        long sum_c_bar = 0;
        long sum_d_bar = 0;

        long rapB = LR;
        long[] raplb = new long[w + 1 - v];
        long[] rapub = new long[w + 1 - v];
        List<Function> rapObj = obj.subList(v, w + 1);

        for (int i = v; i < w + 1; i++) {

            raplb[i - v] = lbVar[i]; //Math.max(lbCopyMDA[i], Math.min(lbVar[i], ubCopyMDA[i]));
            rapub[i - v] = ubVar[i]; //Math.min(ubCopyMDA[i], Math.max(ubVar[i], lbCopyMDA[i]));

        }


        //case three: the RAP problem is feasible under the original bound
        RAP rap = new RAP(rapObj, rapB, raplb, rapub);

        ResultTypeRAP res = rap.solveRAP();
        //ResultTypeRAP res = rap.solveRAPLinear();
        // It is possible that the problem is infeasible.
        if (!res.feasible) {
            //System.out.println("Subproblem" + v + " " + w + "Infeasible");
        }
        return res.sol;
    }

}   

