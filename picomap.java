import java.util.*;
import java.lang.*;
import java.io.*;

public class picomap {
	public static void main(String[] args){ 
		String fm_file = args[0];
		String read_file = args[1];
		int mismatch_penalty = Integer.parseInt(args[2]);
		int gap_penalty = Integer.parseInt(args[3]);
		String output = args[4];

		//read in fm index
		int[] aTally = new int[1];
		int[] cTally = new int[1];
		int[] gTally = new int[1];
		int[] tTally = new int[1];
		//ArrayList<Character> bwtIndex = new ArrayList<Character>();
		StringBuffer genome = new StringBuffer("");
		int[] suffixes = new int[1];


		try {
			// Open the binary file for reading
			FileInputStream fileInputStream = new FileInputStream(fm_file);
			DataInputStream dataInputStream = new DataInputStream(fileInputStream);
			try {
				int stringLen = dataInputStream.readInt();
				for(int i = 0; i < stringLen; i++) {
					genome.append(dataInputStream.readChar());
				}

				int suffixSize = dataInputStream.readInt();
				suffixes = new int[suffixSize];
				for(int i = 0; i < suffixSize; i++) {
					suffixes[i] = dataInputStream.readInt();
				}
				aTally = new int[suffixSize];
				cTally = new int[suffixSize];
				gTally = new int[suffixSize];
				tTally = new int[suffixSize];

				int bwtSize = dataInputStream.readInt();
				bwtIndex = new char[bwtSize];
				for(int i = 0; i < bwtSize; i++) {
					bwtIndex[i] = dataInputStream.readChar();
					//bwtIndex.add(dataInputStream.readChar());
				}

				//readtallies
				for(int i = 0; i < 4; i++) {
					for(int j = 0; j < suffixes.length; j++) {
						if(i == 0) {
							aTally[j] = dataInputStream.readInt();

						} else if(i == 1) {
							cTally[j] = dataInputStream.readInt();

						} else if(i == 2) {
							gTally[j] = dataInputStream.readInt();

						} else if(i == 3) {
							tTally[j] = dataInputStream.readInt();

						}
					}
				}

			} catch (IOException e) {
			}
			dataInputStream.close();
		} catch (IOException e) {
		}


		try(FileReader readReader = new FileReader(read_file)) {
			try (BufferedReader br = new BufferedReader(readReader);
					FileWriter fileWriter = new FileWriter(output);
					BufferedWriter bufferedWriter = new BufferedWriter(fileWriter);
					PrintWriter printWriter = new PrintWriter(bufferedWriter)) {
				String line = "";
				String read_name= "";
				String read = "";
				String reference = "";
				while((line = br.readLine()) != null){		
					if(line.startsWith(">")) {
						read_name = line.substring(1);
						read = br.readLine();
						reference = br.readLine();
					}

				}
				bufferedWriter.close();
				fileWriter.close();
				br.close();
				problemRead.close();
			} catch(IOException e) {

			}
		} catch(IOException e) {

		}
	}
}

for(String name: queryList.keySet()) { 
	int match_len = 0;
	String query = queryList.get(name);
	int top = 1;
	int bottom = bwtIndex.size() - 1;
	for(int i = query.length() - 1; i >= 0; i--) {
		char curr  = query.charAt(i);
		if(curr == 'A') {
			// gets rank of curr char
			top = aTally[top];
			bottom = aTally[bottom];

		} else if(curr == 'C') {
			top = cTally[top] + aTotal;
			bottom = cTally[bottom] + aTotal;

		} else if(curr == 'G') {
			top = gTally[top] + aTotal + cTotal;
			bottom = gTally[bottom] + aTotal + cTotal;

		} else if(curr == 'T') {
			top = tTally[top] + aTotal + cTotal + gTotal;
			bottom = tTally[bottom] + aTotal + cTotal + gTotal;

		}

		if(top == bottom && i != 0) {
			char next = query.charAt(i - 1);
		
				if(bwtIndex.get(i - 1) != next) {
					match_len = 0;
					break;
				}
			
		}
		match_len++;
	}
}

					// construct 2d matrix
					int[][] matrix = new int[read.length() + 1][reference.length() + 1];
				
					for(int j = 0; j <= reference.length(); j++) {
						matrix[0][j] = -1*gap_penalty * j;
					}
					for(int i = 0; i <= read.length(); i++){ 
						matrix[i][0] = -1*gap_penalty * i;
					}
					// fill 2d matrix
					for(int j = 1; j <= reference.length(); j++) {
						for(int i = 1; i <= read.length(); i++) {
							int cost = 0;
							if(read.charAt(i - 1) != reference.charAt(j - 1)) {
								cost = mismatch_penalty;
							} 
							matrix[i][j] = Math.max(Math.max(matrix[i - 1][j] - gap_penalty, matrix[i][j - 1] - gap_penalty), matrix[i - 1][j - 1] - cost);
						}
					}
					// backtrace
					StringBuffer cigar = new StringBuffer("");
					Object[] cig = new Object[read.length() * 2];
					int cigIndex = 2;
					cig[0] = 3;
					cig[1] = "L";
					int i = read.length();
					int j = reference.length();
					int maxj = j;
					// find starting point to allow gaps after reference
					int score = matrix[i][maxj];
					//global alignment condition must reach very corner
					while(i > 0 || j > 0) {
						// fitting alignment condition allows gaps in read
						// resulting i is starting point
						if(i > 0 && matrix[i - 1][j] - gap_penalty == matrix[i][j]) {
							if(cig[cigIndex -1].equals("I")) {
								cig[cigIndex - 2] = (int)cig[cigIndex - 2] + 1;
							} else {
								cig[cigIndex] = 1;
								cig[cigIndex + 1] = "I";
								cigIndex += 2;
							}
							i--;
							//add insertion

						} else if(j > 0 && matrix[i][j - 1] - gap_penalty == matrix[i][j]) {
							if(cig[cigIndex -1].equals("D")) {
								cig[cigIndex - 2] = (int)cig[cigIndex - 2] + 1;
							} else {
								cig[cigIndex] = 1;
								cig[cigIndex + 1] = "D";
								cigIndex += 2;
							}
							j--;
							// add deletion

						} else if(i > 0 && j > 0 && (matrix[i - 1][j - 1] - mismatch_penalty ==  matrix[i][j] || matrix[i - 1][j - 1] == matrix[i][j])) {
							if(matrix[i - 1][j - 1] == matrix[i][j]) {
								//match
								if(cig[cigIndex - 1].equals("=")) {
									cig[cigIndex - 2] = (int)cig[cigIndex - 2] + 1;
								} else {
									cig[cigIndex] = 1;
									cig[cigIndex + 1] = "=";
									cigIndex += 2;
								}

							} else {
								//mismatch
								if(cig[cigIndex -1].equals("X")) {
									cig[cigIndex - 2] = (int)cig[cigIndex - 2] + 1;
								} else {
									cig[cigIndex] = 1;
									cig[cigIndex + 1] = "X";
									cigIndex += 2;
								}
							}
							i--;
							j--;

						} else if(true) {
							System.out.println(i);
							System.out.println(j);
						}
					}
					for(int index = cigIndex - 1; index > 2; index -= 2) {
						cigar.append(String.valueOf(cig[index - 1]));
						cigar.append(String.valueOf(cig[index]));
					}
				}
				