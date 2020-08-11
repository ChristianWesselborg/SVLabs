package pdbfiler;
import pdbfiler.*;

public class PDBInput{

	public void inputFile(String testfile){
		System.out.println(testfile);
	}

	public void run(String passfile){
		Writer go = new Writer();
		go.Writseq(passfile);
	}

	public static void main(String args[]){
		try{
			StringInput input = new StringInput();
			String file = input.userInput();
			
			System.out.println("Entered file: " + file);
			System.out.println("Performing initial checks now...");

			PDBInput pdbinput = new PDBInput();
			pdbinput.inputFile(file);

			FileExists exists = new FileExists();
			if(exists.exist(file) == true){
				System.out.println("Success: file checked");
			}
			else{
				System.out.println("Get file input failed");
				throw new RuntimeException();
			}
			
			PDBCheck check = new PDBCheck();
			if(check.ext(file) == true){
				System.out.println("Success: .pdb file recognized");
			}
			else{
				System.out.println(".pdb check failed");
				throw new RuntimeException(".pdb check failed");
			}

			PDBInput start = new PDBInput();
			start.run(file);
		}
		catch(Exception e){
			System.out.println("Program failed");
		}
	}
	
}
