package jmetal.init;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.StringReader;
import java.nio.charset.Charset;

import javax.swing.JFileChooser;


public class ExtractResults {

	static int numberOfFiles = 35;
	static int numberOfLinesForInformation = 58;
	static int numberOfGenerations= 100;
	static double numberOfEvaluations= 10000;
	
	static double hypervolume=0;
	static int numberOfExactEvaluations = 0;
	static double numberOfSolutionsOnTheFront = 0;
	static int executionTimeHours = 0;
	static int executionTimeMinutes = 0;
	static int executionTimeSeconds = 0;
	static int executionTimeTransformedToSeconds = 0;
	
	public static void main(String[] args) throws IOException {

		String wd = System.getProperty("user.dir");
		wd = "D:\\Doktorat\\Testiranje\\New old testing\\Rezultati";
		
		for(int i=0; i<numberOfFiles; i++)
		{		
			JFileChooser fc = new JFileChooser(wd);
			int rc = fc.showDialog(null, "Select Data File");
			
			if (rc == JFileChooser.APPROVE_OPTION)
			{
				File file = fc.getSelectedFile();
				//getResultFromTheGPDEMOFile(file);
				//getResultFromTheDEMOFile(file);
				getResultFromTheGECFile(file);
				wd = file.getParent();
			}
			else
			{
				numberOfFiles = i;
				System.out.println("File chooser cancel button clicked");
			}
		}	

		executionTimeTransformedToSeconds = 3600 * executionTimeHours + 60*executionTimeMinutes + executionTimeSeconds;
		executionTimeTransformedToSeconds = executionTimeTransformedToSeconds/numberOfFiles;
		
		executionTimeHours = executionTimeTransformedToSeconds / 3600;
		executionTimeMinutes = (executionTimeTransformedToSeconds % 3600)/60;
		executionTimeSeconds = (executionTimeTransformedToSeconds % 3600)%60;
		
		System.out.println("Exact evaluations: " + numberOfExactEvaluations/numberOfFiles);
		System.out.println("Execution time: " + executionTimeHours + ":" + executionTimeMinutes + ":" + executionTimeSeconds);
		System.out.println("Hypervolume: " + hypervolume/numberOfFiles);
		System.out.println("Number of solutions on the front: " + numberOfSolutionsOnTheFront/numberOfFiles);
		System.out.println("Average execution of solution: " + executionTimeTransformedToSeconds/numberOfEvaluations);
	}

	private static void getResultFromTheGPDEMOFile(File file) throws IOException {
	  
		InputStream fis = new FileInputStream(file.getAbsolutePath());
		BufferedReader br  = new BufferedReader(new InputStreamReader(fis, Charset.forName("UTF-8")));
		int numberOfAdditionalEvaluations = 0;
		int lineNumber = 0;
		String line;
		String[] words;
		String temp;
		String time;
	  
		for(lineNumber=0; lineNumber<numberOfLinesForInformation + numberOfGenerations; lineNumber++)
			br.readLine();
		
		//last line of results
		line = br.readLine();
		words = line.split("\t");
		//hypervolume = hypervolume + Double.valueOf(words[6]);
		numberOfExactEvaluations = numberOfExactEvaluations + Integer.valueOf(words[11]);
		
		line = br.readLine();
		temp = line.substring(line.indexOf(':') + 2);
		numberOfAdditionalEvaluations = numberOfAdditionalEvaluations + Integer.valueOf(temp);
		
		numberOfExactEvaluations = numberOfExactEvaluations + numberOfAdditionalEvaluations;
		
		line = br.readLine();
		hypervolume = hypervolume + Double.valueOf(line.substring(line.indexOf(':') + 2));
		line = br.readLine();
		line = br.readLine();
		time = line.substring(line.indexOf(':') + 2);
		executionTimeHours = executionTimeHours + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeMinutes = executionTimeMinutes + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeSeconds = executionTimeSeconds + Integer.valueOf(time);
		
		
  }


	private static void getResultFromTheDEMOFile(File file) throws IOException {
	  
		InputStream fis = new FileInputStream(file.getAbsolutePath());
		BufferedReader br  = new BufferedReader(new InputStreamReader(fis, Charset.forName("UTF-8")));
		int numberOfAdditionalEvaluations = 0;
		int lineNumber = 0;
		String line;
		String[] words;
		String temp;
		String time;
	  
		for(lineNumber=0; lineNumber<numberOfLinesForInformation + numberOfGenerations; lineNumber++)
			br.readLine();
		
		//last line of results
		line = br.readLine();
		words = line.split("\t");
		//hypervolume = hypervolume + Double.valueOf(words[6]);
		hypervolume = hypervolume + Double.valueOf(words[6]);

		line = br.readLine();
		line = br.readLine();
		time = line.substring(line.indexOf(':') + 2);
		executionTimeHours = executionTimeHours + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeMinutes = executionTimeMinutes + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeSeconds = executionTimeSeconds + Integer.valueOf(time);
		
		
  }


	private static void getResultFromTheGECFile(File file) throws IOException {
	  
		InputStream fis = new FileInputStream(file.getAbsolutePath());
		BufferedReader br  = new BufferedReader(new InputStreamReader(fis, Charset.forName("UTF-8")));
		int numberOfAdditionalEvaluations = 0;
		int lineNumber = 0;
		String line;
		String[] words;
		String temp;
		String time;
	  
		for(lineNumber=0; lineNumber<numberOfLinesForInformation + numberOfGenerations; lineNumber++)
			br.readLine();
		
		//last line of results
		line = br.readLine();
		line = br.readLine();	//additional for recalculating hypervolume
		words = line.split("\t");
		//hypervolume = hypervolume + Double.valueOf(words[6]);
		hypervolume = hypervolume + Double.valueOf(words[6]);
		numberOfSolutionsOnTheFront = numberOfSolutionsOnTheFront + Double.valueOf(words[5]);
		numberOfExactEvaluations = numberOfExactEvaluations + Integer.valueOf(words[11]);;
		
		line = br.readLine();
		line = br.readLine();
		time = line.substring(line.indexOf(':') + 2);
		executionTimeHours = executionTimeHours + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeMinutes = executionTimeMinutes + Integer.valueOf(time.substring(0, time.indexOf(':')));
		
		time = time.substring(time.indexOf(':') + 1);
		executionTimeSeconds = executionTimeSeconds + Integer.valueOf(time);
		
		
  }

}
