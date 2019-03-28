import java.io.BufferedWriter;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Scanner;

public class ourSatellite
{
	static double Pi, c, R, S;
	static double[][] V;
	static double[] Xv, Xs;
	static double Tv, Ts;
	static String path;
	static ArrayList<String> output;
	
	public static void main(String[] args) throws IOException
	{
		
		Ts = 0.0d;
		Xs = null;
		Initialize();
		ReadInput();
		WriteOutput();
	}
	
	private static void ReadInput() throws IOException
	{
		Scanner In = new Scanner(System.in);
		String s;
		while((s = In.nextLine()) != null) { ProcessLine(s); }
		In.close();
	}

	private static void ProcessLine(String s)
	{
		String[] split = s.split(" ");
		double[] pos = PositionInCartesian(split);
		Tv = Double.parseDouble(split[0]);
		Xv = pos;
		// Gets the satellites that are visible for the position
		boolean[] vis = VisibleCheck(pos);
		for(int i = 0; i < vis.length; i++)
		{
			if(vis[i])
			{
				Ts = ComputeTs(i);
				Xs = Xs(i, Ts);
				String ret = i+" "+Ts+" "+Xs[0]+" "+Xs[1]+" "+Xs[2];
				System.out.println(ret);
				output.add(ret);
				if(Ts>200 && Ts<203)
				{
					System.out.println("Look here!");
					output.add("Look here!");
				}
			}
		}
	}
	
	static void WriteOutput()
	{
		try(Writer writer = new BufferedWriter(new OutputStreamWriter(new FileOutputStream("satellite.log"),"utf-8")))
		{
			for(int i = 0; i<output.size(); i++) { writer.write(output.get(i)+"\n"); }
			writer.close();
		}
		catch(IOException e) { System.out.println("Something went wrong."); }
	}
	
	private static double ComputeTs(int v)
	{
		double t0 = Tv;
		double norm = 0.0d;
		for(int i = 0; i < 3; i++)
		{
			norm += Math.pow(Xs(v,Tv)[i],2);
		}
		t0 = t0 - Math.sqrt(norm)/c;
		
		return NewtonsMethod(v,t0,0);
	}
	
	private static void Initialize()throws IOException
	{
		V = new double[24][9];
		output = new ArrayList<String>();
		Scanner scan = new Scanner(new File("data.dat")).useDelimiter("\n");
		String[] l = new String[220];
		int counter = 0;
		while(scan.hasNext())
		{	
			l[counter] = scan.next();;
			counter++;
		}
		scan.close();
		Pi = Double.parseDouble(l[0].substring(0, 27));
		c = Double.parseDouble(l[1].substring(0, 27));
		R = Double.parseDouble(l[2].substring(0, 27));
		S = Double.parseDouble(l[3].substring(0, 27));
		
		counter = 4;
		for(int v = 0; v < 24; v++)
		{
			for(int i = 0; i < 9; i++)
			{
				V[v][i] = Double.parseDouble(l[counter].substring(0, 27));
				counter++;
			}
		}
		scan.close();
	}
	
	// Computations	
	private static double[] Xs(int v, double t)
	{
		double x = (R+V[v][7])*(V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double y = (R+V[v][7])*(V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		double z = (R+V[v][7])*(V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]) + V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	private static double[] Xs1(int v, double t)
	{
		double x = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][3]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][0]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double y = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][4]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][1]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		double z = 2*Pi*(1/V[v][6])*(R+V[v][7])*(V[v][5]*Math.sin((2*Pi*t)/V[v][6]+V[v][8])-V[v][2]*Math.cos((2*Pi*t)/V[v][6]+V[v][8]));
		return new double[]{x,y,z,t};
	}
	private static double f(int v, double t)
	{
		double calc = -1*c*c*(Tv-t)*(Tv-t);
		for(int i = 0; i<3; i++)
		{
			calc += (Xs(v,t)[i]-Xv[i])*(Xs(v,t)[i]-Xv[i]);
		}
		return calc;
	}
	private static double f1(int v, double t)
	{
		double retval = 2*c*c*(Tv-t);
		for(int i = 0; i<3; i++)
		{
			retval += 2*(Xs(v,t)[i]-Xv[i])*Xs1(v,t)[i];
		}
		return retval;
	}
	private static double NewtonsMethod(int v, double Tk, int depth)
	{
		double iter = Tk;
		iter = Tk-(f(v,Tk)/f1(v,Tk));
		if(iter-Tk < 0.01/c) { return iter; }
		else if(depth>=9) { return -1; }
		else { return NewtonsMethod(v, iter, depth+1); }
	}

	// Convert position to Cartesian coordinates
	private static double[] PositionInCartesian(String[] a)
	{
		double 	t 		= 	Double.parseDouble(a[0]);
		double 	theta 	= 	Integer.parseInt(a[4])*DegreesToRadians(new String[]{a[1],a[2],a[3]});
		double 	phi		= 	Integer.parseInt(a[8])*DegreesToRadians(new String[]{a[5],a[6],a[7]});
		double 	h 		= 	Double.parseDouble(a[9]);
		double 	x 		= 	(R + h) * Math.cos(theta) * Math.cos(phi);
		double 	y 		= 	(R + h) * Math.cos(theta) * Math.sin(phi);
		double 	z 		= 	(R + h) * Math.sin(theta);
		double  alpha   =   (2*Pi*t)/S;
		
		return new double[]{Math.cos(alpha)*x-Math.sin(alpha)*y,Math.sin(alpha)*x+Math.cos(alpha)*y,z,t};
	}
	// Convert degrees to Radians
	private static double DegreesToRadians(String[] a)
	{
		return 2*Pi*(Integer.parseInt(a[0])/360.0d + Integer.parseInt(a[1])/(360.0d*60) + Double.parseDouble(a[2])/(360*60*60));
	}
	
	// Checks to see if satellite visible
	private static boolean[] VisibleCheck(double[] a)
	{
		boolean[] vis = new boolean[24];
		double[] satPos;
		for(int i = 0; i < 24; i++)
		{
			satPos = Xs(i, a[3]);
			vis[i] = 2*a[0]*(satPos[0]-a[0]) + 2*a[1]*(satPos[1]-a[1]) + 2*a[2]*(satPos[2]-a[2]) > 0;
		}
		return vis;
	}
}
