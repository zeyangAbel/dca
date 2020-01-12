/******************************************************************************
* This is a class for storing the solution to the subproblems of DRAP-NC in the MDA method
****************************************************************************/

public class ResultTypeMDA {
	long[] aa;
	long[] ab;
	long[] ba;
	long[] bb;

	public ResultTypeMDA(long[] aa, long[] ab, long[] ba, long[] bb) {
		this.aa = aa;
		this.ab = ab;
		this.ba = ba;
		this.bb = bb;
	}
} 