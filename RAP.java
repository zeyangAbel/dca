/******************************************************************************
 *  Compilation:  javac RAP.java
 *  Execution:    java RAP
 *  
 *	@author Zeyang Wu University of Minnesota
 *
 * This is a class for discrete simple resource allocation (DRAP) with separable convex objectives. 
 *      
 * Method: solveRAP()
 *         
 * The scaling algorithm in Hochbaum 1994 paper is implemented. The worse cases running time is 
 * O(nlognlonB) where n is the number of periods and B is total amount of demands.  
 * subroutine: greedyRAP(int, long[])
 *		  
 * Method: solveRAPLinear()
 * An O(n) time algorithm for DRAP with linear objective. The exact procedure is given by Hochbaum in 1995:
 * 		About strongly polynomial time algorithms for quadratic optimization over sub-modular constraints.
 * 
*/

import java.util.*;

public class RAP {
	//list of function oracles: The user should write their own function classes which extends the abstract class Function .
	List<Function> obj;
	long B;
	long[] lbVar;
	long[] ubVar;
	int dimension;
	long scaleFactor;

	public RAP(List<Function> obj, long B, long[] lbVar, long[] ubVar) {
		this.obj = obj;
		this.B = B;
		this.lbVar = lbVar;
		this.ubVar = ubVar;
		this.dimension = lbVar.length;
		this.scaleFactor = 1;		
	}


	/**
     * Solve operations.
     * With this operation you can solve the simple resource allocation problem with separable convex objectives
     * The is an implementation of Hochbaum 1994 Lower and Upper Bounds for the Allocation Problem and Other Nonlinear Optimization Problems
     * <p>
     * Time-Complexity: O(n log n log(B)) 
     * The running time should be improved to O(n log(B)) if the CUT procedure is implemented in Solve()
     * However, the implementation is non-trivial.
     * @param no param
     */
	public ResultTypeRAP solveRAP() {
		if (dimension == 1) {
			if (B <= ubVar[0] && B >= lbVar[0]) {
				long[] x = new long[] {B};
				return new ResultTypeRAP(true, x);
			} else {
				return new ResultTypeRAP(false, null);
			}
		}
		long[] x = new long[dimension];
		//create a feasible solution.
		for (int i = 0; i < dimension; i++) {
				x[i] = lbVar[i];
		}

		//Step size s
		long s = (long) Math.ceil(((double) B) / dimension / 2);
		while (s > 1) {
			//greedy step with size s
			greedyRAP(s, x, B);

			//undo the last step of greedy(s)
			for (int i = 0; i < dimension; i++) {
				//avoid infeasibility				
				if (x[i] - s < lbVar[i]) {
					x[i] = lbVar[i];
					continue;
				}				
				x[i] = x[i] - s;
			}
 
 			//New step size
			if (s % 2 == 0) {
				s = s / 2;
			} else {
				s = s / 2 + 1;
			}
		}
        
        //return
        return greedyRAP(1, x, B);		
	}

	//The DataType class represents a resource with its id and unit allocation cost
    private class DataType{
    	int id;
    	double value;
    	public DataType(int id, double value) {
    		this.id = id;
    		this.value = value;
    	}
    }

    private class DataTypeComparator implements Comparator<DataType> {
    	public int compare(DataType a, DataType b) {    		
    		if (a.id == b.id) {
    			return 0;
    		}    		
    		if (a.value > b.value) {
    			return 1;
    		}
    		return -1;
    	}
    }


    /**
     * Greedy Algorithm(Marginal Allocation Algorithm) with step size s.
     * Balanced Binary Search tree is used to store the activities id and unit allocation cost.
     * <p>
     * Time-Complexity: O(m log(n)) where m is the total number of increments made during this process 
     *
     * @param s step size
     * @param x a initial vector can be 0 or any component-wise lower of the optimal solution  
     * @param B total amount of resource
     */

	public ResultTypeRAP greedyRAP(long s, long[] x, long B) {
		boolean[] removedIndex = new boolean[dimension];
		long numOfRemovedIndex = 0;
		for (int i = 0; i < dimension; i++) {
			B -= x[i];
		}

		TreeSet<DataType> heap = new TreeSet<>(new DataTypeComparator());
		DataType[] indexList = new DataType[dimension]; 
		for (int i = 0; i < dimension; i++) {
			double increase =  obj.get(i).getValue(((double) (x[i] + 1))/scaleFactor) - obj.get(i).getValue( ((double) x[i]) / scaleFactor);
			DataType temp = new DataType(i, increase);
			heap.add(temp);
			indexList[i] = temp;
		}

		while (B >= 1 && numOfRemovedIndex < dimension) {
			//find the minimum increase
			int minIndex = -1;
			double minIncrease = Long.MAX_VALUE;

			//O(log n) operation, which can be improved to O(1) amortized time by CUT
			//The CUT procedure Not being implemented yet
			minIndex = heap.pollFirst().id;


			//increase x[minIndex]
			//feasibility check which takes O(1) time
			if (x[minIndex] + 1 > ubVar[minIndex]) {
				removedIndex[minIndex] = true;
				numOfRemovedIndex++;
			} else if (s > 1 && (x[minIndex] + s > ubVar[minIndex] || B < s)) {
				removedIndex[minIndex] = true;
				numOfRemovedIndex++;
				//Errata in the 2008 paper 
				long increaseUnit = Math.min(ubVar[minIndex] - x[minIndex], B);
				x[minIndex] += increaseUnit;
				B -= increaseUnit;
			} else {
				x[minIndex] += s;
				B -= s;
				indexList[minIndex].value = obj.get(minIndex).getValue(((double) (x[minIndex] + 1))/scaleFactor) - obj.get(minIndex).getValue( ((double) x[minIndex]) / scaleFactor);
				heap.add(indexList[minIndex]);
			}
		}
		
		//System.out.println(B);
		if (s > 1 || B == 0) {
			return new ResultTypeRAP(true, x);
		} else {
			return new ResultTypeRAP(false, null);
		}
	}

	/**
     * SolveRAPLinear().
     * With this operation, we can solve the RAP problem in O(n) time when the objective function is linear.
     * <p>
     * This algorithm is first derived by Brucker in 1984 to solve a linear time algorithm for the continuous convex quadratic Knapsack problem
     * The problem studied is a general version of simple resource allocation problem. 
     * The exact procedure is given by Hochbaum in 1995:
     * About strongly polynomial time algorithms for quadratic optimization over sub-modular constraints 
     * Time-Complexity: O(n) 
     * 
	 */

	public ResultTypeRAP solveRAPLinear() {
		for (int i = 0; i < dimension; i++) {
			B -= lbVar[i];
			ubVar[i] -= lbVar[i];
		}
		ResultTypeRAP res = solveRAPLinear(null);
		for (int i = 0; i < dimension; i++) {
			res.sol[i] += lbVar[i];
		}
		return res;
	}
	 
	public ResultTypeRAP solveRAPLinear(String[] args) {
	 	//break points
	 	DataType[] breakPoint = new DataType[dimension];
	 	for (int i = 0; i < dimension; i++) {
	 		//For linear functions, we can get its coefficient by getValue(1);
	 		breakPoint[i] = new DataType(i, obj.get(i).getValue(2) - obj.get(i).getValue(1));
	 	}

	 	ResultTypeRAP res = new ResultTypeRAP(true, new long[dimension]);

	 	//bisection to find dual break point;
	 	int start = 0;
	 	int end = dimension - 1;
	 	long sum = 0;
	 	while(start + 1 < end) {
	 		int mid = start + (end - start) / 2;

	 		this.kthSmallest(breakPoint, start, end, mid - start + 1);
	 		
	 		//prefix sum to reduce time complexity from O(n^2) to O(n)
	 		for (int i = start; i < mid + 1; i++) {
	 			int index = breakPoint[i].id;
	 			sum += ubVar[index];
	 		}

	 		if (sum > B) {
	 			end = mid;
	 		} else {
	 			start = mid;
	 		}

	 		//reset prefix sum so that it is not affected by the function kthsmallest().
	 		if (start == 0) {
	 			sum = 0;
	 		} else {
	 			for (int i = mid; i >= start; i--) {
	 				int index = breakPoint[i].id;
	 				sum -= ubVar[index];
	 			}
	 		}

	 	}

	 	sum = 0;
	 	for (int i = 0; i < start; i++) {
	 		int index = breakPoint[i].id;
	 		res.sol[index] = ubVar[index];
	 		sum += ubVar[index];	 			
	 	}

	 	//trivial cases (two break points)
	 	//In case that The bisection is not performed for (end - start + 1) = 2.
	 	//This is essential
	 	if (breakPoint[start].value > breakPoint[end].value) {
	 		swap(breakPoint, start, end);
	 	}
	 	//Trivial cases
 		int index = breakPoint[start].id;
	 	if (sum + ubVar[index] > B) {
	 		res.sol[index] = B - sum;
	 		sum = B;
		} else {
	 		sum += ubVar[index];
	 		res.sol[index] = ubVar[index];
	 	}
	 	
	 	if (start < end) {
	 		index = breakPoint[end].id;
	 		res.sol[index] = B - sum;
	 	} 	

	 	return res;
	 }


	/**
	* kthSmallest operation: median of the median algorithm
  	* Returns k'th smallest element in arr[l..r] in worst case linear time. 
  	* ASSUMPTION: ALL ELEMENTS IN ARR[] ARE DISTINCT
  	* <p>
  	* Time-Complexity: O(n) where n is the length of the array 
  	* @author adapt from the C++ code for geeksForgeeks
  	* @param array the array 
  	* @param start the start index
  	* @param end the end index
  	* @param k the k-th smallest element in a[s]...a[e]
  	*/
  	public DataType kthSmallest(DataType[] array, int start, int end, int k) {
  		// If k is smaller than number of elements in array
  		if (k <= 0 || k > end - start + 1) {
  			return null;
  		}

  		//number of elements
  		int n = end - start + 1;

  		// Divide arr[] in groups of size 5, calculate median of every group and store it in median[] array.
  		int i = 0;
  		DataType[] median = new DataType[(n+4) / 5];

  		for (i = 0; i < n / 5; i++) {
  			int left = start + 5 * i;
  			int right = left + 4;

  			Arrays.sort(array, left, right + 1, new DataTypeComparator());
  			int mid = left + (right - left) / 2;
      		median[i] = array[mid];
      	}

      	//For the last few elements
        if (i * 5 < n) {
      		Arrays.sort(array, start + 5*i, end + 1, new DataTypeComparator());
      		median[i] = array[start + 5*i + n%5/2];
            i++;
      	}

      	// Find median of all medians using recursive call.
    	// If median[] has only one element, then no need
    	// of recursive call
    	DataType medOfMed = (i == 1) ? median[i - 1] : kthSmallest(median, 0, i - 1, i / 2);

    	int pos = partition(array, start, end, medOfMed);

    	// If position is same as k
    	if (pos - start == k-1) {
      		return array[pos];
    	}

    	// If position is more, recur for left
    	if (pos - start > k - 1) {
      		return kthSmallest(array, start, pos - 1, k);
    	}
    
    	// Else recur for right subarray
    	return kthSmallest(array, pos + 1, end, k - pos + start -1);
  	}


  	private static void swap(DataType[] arr, int i, int r) {
  		DataType temp = arr[i];
  		arr[i] = arr[r];
  		arr[r] = temp;
  	}


  	// It searches for x in arr[l..r], and partitions the array 
	// around x.
	private static int partition(DataType[] arr, int l, int r, DataType x) {
    	// Search for x in arr[l..r] and move it to end
    	int i;
    	for (i = l; i < r; i++) {
        	if (arr[i].value == x.value) {
        		break;
        	}
        }

    	swap(arr, i, r);
 
    	// Standard partition algorithm
    	i = l;
    	for (int j = l; j <= r - 1; j++) {
        	if (arr[j].value <= x.value) {
        	    swap(arr, i, j);
            	i++;
        	}
    	}
    	swap(arr, i, r);

    	return i;
	}

}