package jmetal.init;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.DataInputStream;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStreamReader;


public class EKGWraper {

	/**
	 * @param args
	 * @throws IOException 
	 */
	public static void main(String[] args) throws IOException {
		
	/*	if(args.length != 8)
		{
			System.out.println("Write 8 input variables");
			System.exit(0);
		}*/
		
		String simulatorCommand = System.getProperty("user.dir") + "\\singleRun.bat";
		File simmulatorDirectory = new File(System.getProperty("user.dir")); // + "\\EKGSimulator\\");	//"D:\\IJS\\EKG_simulator_MDepolli\\");
		File outputFile = new File(System.getProperty("user.dir") + "\\izhod.txt");
		File simOut = new File(System.getProperty("user.dir") + "\\simOut.txt");
		File simIn = new File(System.getProperty("user.dir") + "\\simIn.txt");
		String objectives = "";
		String violation = "";
		String[] variables = new String[8];
		
		try
		{
			FileInputStream fstream = new FileInputStream(simIn);
		  DataInputStream in = new DataInputStream(fstream);
		  BufferedReader br = new BufferedReader(new InputStreamReader(in));
 
		  for(int i=0; i<8;i++)
		  	variables[i] = br.readLine();
		}		
		catch (Exception e){//Catch exception if any
			System.out.println("Write 8 input variables");
		  System.err.println("Error: " + e.getMessage());
		}
		
		try{
		  // Create file 
		  FileWriter fstream = new FileWriter(simulatorCommand);
		  BufferedWriter out = new BufferedWriter(fstream);
		  out.write("ekgsim.exe -sim ");
		  out.write(variables[0]);
		  for(int i=1; i< variables.length; i++)
		  	out.write("," + variables[i]);
		  //Close the output stream
		  out.close();
		  fstream.close();
		}
		catch (Exception e){//Catch exception if any
		  System.err.println("Error: " + e.getMessage());
		}
		
		try
		{
			ProcessBuilder pp = new ProcessBuilder(simulatorCommand);
    	pp.directory(simmulatorDirectory);
    	pp.redirectOutput(outputFile);
    	Process process = pp.start();
    	process.waitFor();			
		} 
		catch (IOException | InterruptedException e) {	   
	    e.printStackTrace();
    }
		finally	{		}
		
		try{
			  FileInputStream fstream = new FileInputStream(outputFile);
			  DataInputStream in = new DataInputStream(fstream);
			  BufferedReader br = new BufferedReader(new InputStreamReader(in));
			  String strLine;
	
			  while ((strLine = br.readLine()) != null)   {
			  	if(strLine.contains("criteria ="))
			  	{
			  		objectives = strLine.substring(strLine.indexOf('<')+1, strLine.indexOf('>'));
			  		violation = strLine.substring(strLine.indexOf("violation") + 12);
			  		break;
			  	}
	
			  }
			  //Close the input stream
			  in.close();
			  fstream.close();
		}			
		catch (Exception e){//Catch exception if any
			  System.err.println("Error: " + e.getMessage());
		}
		
		String[] tableOfObjectives = objectives.split(",");
		FileWriter fstream = new FileWriter(simOut);
	  BufferedWriter out = new BufferedWriter(fstream);
	  
	  out.write("1");
	  out.newLine();
	  double distanceToFisibility =  Double.parseDouble(violation);
		if(distanceToFisibility != 0)
			out.write(" 1");
		else
			out.write(" 0");
		
		out.newLine();
		
		for(int i=0; i<tableOfObjectives.length; i++)
		{
			out.write(" " + tableOfObjectives[i]);
			out.newLine();
		}
	
		out.close();
	  fstream.close();
	}
}
