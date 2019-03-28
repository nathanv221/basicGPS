//package gps_package;

import java.util.Scanner;
import java.util.ArrayList;


public class ourReceiver
	{
	static double[][] data; //only changed in GetInput() but used in almost
				//everything, and was getting annoying to pass around
				//holds all satellite data
	//static int satnum;		//only changed in GetInput() but used everywhere
					//keeps track of how many satellites we are getting data from
	static double c = Double.parseDouble("2.99792458E+08");
	static double R = 6367444.5;
							//These never change
	static int chunkSize;
	static int location = 0;
	static boolean debug =false; //just for debuging 
	
	
	/*
	 * This program will take the satellite arguments and return the vehicle arguments
	 * the main method will be used to call each method to avoid nesting where possible
	 * and have the final print command
	 */
	public static void main(String[] args) 
	{
		//int i = 0;
		//allChunkSize = 0;
		ArrayList<Double> input = GetInput();
		int dataSize = input.size();
		
		while (location < input.size())
		{
			UseInput(input);
			GetLocation();
		}
	}
	public static void GetLocation()
	{
		double[] xv = Newtons(new double[] {0,0,0}, 0);
		double tv = TimeAt(xv);
		//tv = ((int) ((tv*100.0) + ((tv <0.0) ? -0.5:0.5)))/100.0;
		String pos = CartToGlobe(xv[0], xv[1], xv[2], tv);
		System.out.println(tv + " " + pos);
		
	}
	static double TimeAt(double[] xv)
	{
		return data[0][0] + (distance(xv[0], data[0][1], 
						xv[1], data[0][2], xv[2], data [0][3]))/c;
	}
	
	/*
	 * a simple scanner that puts the input into an arraylist
	 * using arraylist because there could be more than 4 satellites
	 * then puts the input into a double array for storage
	 */
	public static ArrayList<Double> GetInput() 
	{
		ArrayList<Double> input = new ArrayList<Double>();
		
		//satnum = 0;
		int counter = 0;
		Scanner scan = new Scanner(System.in);
		
		if (debug)
			System.out.println("called get input");

		while (scan.hasNext())
		{
			input.add(scan.nextDouble());
			//counter++; 

			/*
			//when counter=5 we are getting data from a new satellite
			if (counter == 5)
			{
				counter = 0;
				satnum++;
			if (debug)
				System.out.print(satnum + " ");
			}
			*/
		}
		scan.close();
		return input;
		
		
	}
	public static double[][] UseInput(ArrayList<Double> input)
	{
		//now lets put that data somewhere useful


		//this part was a bitch to figure out. checks if the data is related to the previous 
		//data receved this accounts for half our global vars
		chunkSize =1;
		while ((location +1+(5*chunkSize))<input.size() 
			&& Math.abs(input.get(location+1+(5*chunkSize))-input.get(location+1+(5*(chunkSize-1)))) < 0.5)
		{
			chunkSize++;
		}

		//allChunkSize = allChunkSize + chunkSize;
		data = new double[chunkSize][4];
		
		//double currentVeh = input.get(0);
		
		
			for(int j = 0; j < chunkSize; j++)
			{
				// input 0 is the vehicle id
			
				// Input 2 is the time at which the vehicle broadcasts
				data[j][0] = input.get(location + 1);
				// Input 3, 4, 5 are the coordinates (x,y,z) of the satellite
				data[j][1] = input.get(location + 2);
				data[j][2] = input.get(location + 3);
				data[j][3] = input.get(location + 4);
				
				location = location+5;//get next satellites data
			}
		
		return data;
		       
	}


	/*
	 * runs newtons method max of 10 times
	 * returns estimate of location within 1cm
	 * if nonconvergent fails whole program
	 */
	public static double[] Newtons(double[] xn, int count)
	{
		//the single line all this code is based on
		double temp[] = CramersRule(Jacobian(xn), gradf(xn));
		double[] sol = xn.clone();
		for(int i = 0; i<3; i++)
			{ sol[i] = sol[i] - temp[i]; } //newtons method
		if(Math.abs(VectorSubtract(sol,xn)[0])<0.0000001
				&& Math.abs(VectorSubtract(sol,xn)[1])<0.0000001
				&& Math.abs(VectorSubtract(sol,xn)[2])<0.0000001) //off <1cm
			
			{ return sol; } //were done!
		else if(count>10) //if we keep repeating for no reason...
			{ return null; } //just give up
		else
		{
			count++;
			return Newtons(sol, count); //recurse
		}
	}
	/*
	 * runs cramers rule to find the xyz chords
	 */
	static double[] CramersRule(double[][] A, double[] b)
	{
		double coefdet = Determinate3x3(A[0], A[1], A[2]);
		double xdet = Determinate3x3(b, A[1], A[2]);
		double ydet = Determinate3x3(A[0], b, A[2]);
		double zdet = Determinate3x3(A[0], A[1], b);
		
		double x = xdet/coefdet;
		double y = ydet/coefdet;
		double z = zdet/coefdet;
		
		return new double[]{x,y,z};
	}
	/*
	 * simply gets the determinate of a 3x3 matrix
	 */
	public static double Determinate3x3(double[] col1, double[] col2, double[] col3)
	{
		double sol;
		sol = col1[0]*(col2[1]*col3[2] - col2[2]*col3[1])
				-col1[1]*(col2[0]*col3[2] - col2[2]*col3[0])
				+col1[2]*(col2[0]*col3[1] - col2[1]*col3[0]);
							
		return sol;
	}
	
	/*
	 * gets a 3x3 jacobian matrix to use with newtons method
	 */
	public static double[][] Jacobian(double[] x0)
	{
		double[][] jacobian = new double[3][3]; //make a 3x3 matrix
		for(int i = 0; i < 3; i++)
		{
			for(int j = 0; j < 3; j++)
			{
				jacobian[i][j] = JacobianElement(i, j, x0); //this is messy, so it got its own method
						
			}
			
		}
		return jacobian;
		
	}

	private static double JacobianElement(int i, int j, double[] xn) 
	{
		double sum = 0;
		for(int k = 0; k < chunkSize-1; k++) 
		{
			sum = sum + FiniteDiff(k, i, xn) * FiniteDiff(k, j, xn); 
			
		}
		
		return sum*2;
	}
	
	/*
	 * Returns the central finite difference formula approx for 
	 * the Jacobian using the distance from satallite to our newtons
	 * method guess as our h value
	 */
	private static double FiniteDiff(int i, int j, double[] xn)
	{
		double finiteDiff = (xn[j] - data [i][j+1]) 
				/ distance(xn[0], data[i][1], xn[1], data[i][2], xn[2], data[i][3]);
		finiteDiff = finiteDiff - ((xn[j] - data [i+1][j+1]) 
				/ distance(xn[0], data[i+1][1], xn[1], data[i+1][2], xn[2], data[i+1][3]));
		return finiteDiff;
	}

	/*
	 * simple distance formula for a 3 chords
	 */
	private static double distance(double x1, double x2, double y1, double y2, double z1, double z2)
	{
		return Math.sqrt((Square(x1-x2) + Square(y1-y2) + Square(z1-z2)));
	}
	
	/*
	 * javas math square function is ugly, this is easier
	 */
	private static double Square(double x)
	{
		return x*x;
	}
	
	
	static double[] gradf(double[] x)
	{
		//using a double array as an array not a matrix!!!
		double[][] diffs = new double[chunkSize][3];
		
		for(int j = 0; j<chunkSize; j++) 
		{ 
			diffs[j] = VectorSubtract(new double[]{data[j][1],data[j][2],data[j][3]},x);
		}
		
		//create a vector of norms
		double[] normVector = new double[chunkSize];
		
		double norm = 0.0;
		for(int j = 0; j<chunkSize; j++) 
		{ 
			
			for(Double d : diffs[j])
				norm = norm + d*d;
			
			norm = Math.sqrt(norm);
			
			normVector[j]=norm;
		}
		
		
		double[] A = new double[chunkSize-1];
		
		for(int j = 0; j<chunkSize-1; j++)
		{
			A[j] = normVector[j+1]-normVector[j]-c*(data[j][0]-data[j+1][0]);
		}
		
		double[][] XYZ = new double[3][chunkSize-1];
		
		for(int i = 0; i<3; i++)
		{
			for(int j = 0; j<chunkSize-1; j++)
			{
				XYZ[i][j] = diffs[j][i]/normVector[j] - diffs[j+1][i]/normVector[j+1];
			}
		}
		
		double[] sol = new double[3];
		for(int j = 0; j<3; j++)
		{			
			for(int i = 0; i<chunkSize-1; i++) 
			{
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
				sol[j] = A[i]*XYZ[j][i] +sol[j];						///
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////
			}
			sol[j] = sol[j]*2;
		}
		return sol;
	}
	
	/*
	 * super complex code that subtracts 2 vectors
	 */
	static double[] VectorSubtract(double[] a, double[] b)
	{
		if(a.length==b.length)
		{
			double[] sol = new double[a.length];
			
			for(int i = 0; i<a.length; i++)
				sol[i]=a[i]-b[i];
			
			return sol;
		}
		else
			return null;
	}
	public static String CartToGlobe (double x1, double y1, double z, double t)
	{

	// account for the earths rotation
		double rot = -2.0*Math.PI*t/86164.09;
		double x = Math.cos(rot)*x1 - Math.sin(rot) * y1;
		double y = Math.sin(rot)*x1 + Math.cos(rot) * y1;  


		double[]polarchords = CartToPolar(x, y, z);

		double radsLon = polarchords[0];
		double radsLat = polarchords[1];
		double alt = polarchords[2];
		int degreesLat;
		int minutesLat;
		double secondsLat;
		int NS;
		int degreesLon;
		int minutesLon;
		double secondsLon;
		int EW;
		double temp;
		
		//deal with latitude
		
		if (radsLat < 0)
			NS = -1;
		else
			NS = 1;
		
		temp = (Math.abs(radsLat * 180.0/Math.PI));
		degreesLat = (int) Math.floor(temp);
		temp = (temp - (double)degreesLat) * 60.0;
		minutesLat = (int) Math.floor(temp);
		secondsLat = (temp - (double)minutesLat) *60.0;
		
		//deal with longitude
		
		if (radsLon < 0)
			EW = -1;
		else
			EW = 1;
		
		temp = (Math.abs(radsLon * 180.0/Math.PI));
		degreesLon = (int) Math.floor(temp);
		temp = (temp - (double)degreesLon) * 60.0;
		minutesLon = (int) Math.floor(temp);
		secondsLon = (temp - (double)minutesLon) *60.0;
		
		
		return new String(degreesLat + " " + minutesLat + " "
				+ secondsLat + " " + NS + " "+ degreesLon + " "
				+ minutesLon + " " + secondsLon + " " + EW + " " + alt);
		
		//return 0;
	}
	public static double[] CartToPolar(double x, double y, double z)
	{
		double r, azimuthal, polar;
		r = Math.sqrt(Square(x)+Square(y)+Square(z))-R;
		azimuthal = Math.atan2(y,x);
		polar = Math.atan2(z, Math.sqrt(Square(x) + Square(y)));
		
		return new double[] {azimuthal,polar,r};
		
	}
}
