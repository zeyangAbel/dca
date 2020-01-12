/******************************************************************************
* This is a class for storing the solution to an instance of DRAP
****************************************************************************/
public class ResultTypeRAP {
	boolean feasible;
	long[] sol;
	public ResultTypeRAP(boolean feasible, long[] sol) {
		this.feasible = feasible;
		this.sol = sol;
	}
}