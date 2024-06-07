import java.util.*;
import java.io.*;

public class picoindex {
    public static void main(String[] args) {
        String genome = args[0];
		String output = args[1];
		StringBuffer currSequence = new StringBuffer("");
        // read in genome
		try(FileReader genomeRead = new FileReader(genome)) {
			try (BufferedReader br = new BufferedReader(genomeRead)) {
				String line = "";
				while((line = br.readLine()) != null){		
					if(!line.startsWith(">")) {
						currSequence.append(line);
                    }
				}
				br.close();
				genomeRead.close();
			} catch(IOException e) {

			}
		} catch(IOException e) {

		}
		final String sequence = currSequence + "$";
		ArrayList<Integer> suffixes = new ArrayList<Integer>();
		for(int i = 0; i < sequence.length(); i++) {
			suffixes.add(i);
		}
		Collections.sort(suffixes, new Comparator<Integer>() {
            @Override
            public int compare(Integer i, Integer j) {
                while (i < sequence.length() && j < sequence.length()) {
                    if (sequence.charAt(i) != sequence.charAt(j)) {
                        return sequence.charAt(i) - sequence.charAt(j);
                    }
                    i++;
                    j++;
                }
                return i - j;
            }
        });
        // FM Index construction
		// generate bwt index (last column)
		char[] bwtIndex = new char[suffixes.size()];
		for(int i = 0; i < suffixes.size(); i++){ 
			if(suffixes.get(i) == 0){ 
				bwtIndex[i] = '$';
			} else { 
				bwtIndex[i] = sequence.charAt(suffixes.get(i) - 1);
			}
			
		}

		// get tally indexes
		int[] aTally = new int[suffixes.size()];
		int[] tTally = new int[suffixes.size()];
		int[] cTally = new int[suffixes.size()];
		int[] gTally = new int[suffixes.size()];

		if (bwtIndex[0] == 'A') { 
			aTally[0] = 1;

		} else if(bwtIndex[0] == 'T') {
			tTally[0] = 1;
			
		} else if(bwtIndex[0] == 'C') {
			cTally[0] = 1;

		} else if(bwtIndex[0] == 'G') {
			gTally[0] = 1;
		}

		for(int i = 1; i < bwtIndex.length; i++) {
			if(bwtIndex[i] == 'A') {
				aTally[i] = aTally[i - 1] + 1;
				gTally[i] = gTally[i - 1];
				cTally[i] = cTally[i - 1];
				tTally[i] = tTally[i - 1];

			} else if(bwtIndex[i] == 'T') {
				tTally[i] = tTally[i - 1] + 1;
				gTally[i] = gTally[i - 1];
				cTally[i] = cTally[i - 1];
				aTally[i] = aTally[i - 1];  

			} else if(bwtIndex[i] == 'C') {
				cTally[i] = cTally[i - 1] + 1;
				aTally[i] = aTally[i - 1]; 
				gTally[i] = gTally[i - 1]; 
				tTally[i] = tTally[i - 1]; 

			} else if(bwtIndex[i] == 'G') {
				gTally[i] = gTally[i - 1] + 1;
				aTally[i] = aTally[i - 1]; 
				cTally[i] = cTally[i - 1]; 
				tTally[i] = tTally[i - 1]; 

			} else if(bwtIndex[i] == '$') {
				gTally[i] = gTally[i - 1];
				aTally[i] = aTally[i - 1]; 
				cTally[i] = cTally[i - 1]; 
				tTally[i] = tTally[i - 1]; 
			}
		}

		try {
			FileOutputStream fileOutput = new FileOutputStream(output);
			DataOutputStream dataOutput = new DataOutputStream(fileOutput);

			//write genome
			dataOutput.writeInt(sequence.length());
			for(int i = 0; i < sequence.length(); i++){ 
				dataOutput.writeChar(sequence.charAt(i));
			}
			
			//write suffix array
			dataOutput.writeInt(suffixes.size());
			for(int i = 0; i < suffixes.size(); i++){ 
				dataOutput.writeInt(suffixes.get(i));
			}
			// write bwtIndex
			dataOutput.writeInt(bwtIndex.length);
			for(int i = 0; i < bwtIndex.length; i++){ 
				dataOutput.writeChar(bwtIndex[i]);
			}
			
			//write tallies
			for(int i = 0; i < 4; i++) {
				for(int j = 0; j < aTally.length; j++) {
					if(i == 0) {
						dataOutput.writeInt(aTally[j]);

					} else if(i == 1) {
						dataOutput.writeInt(cTally[j]);

					} else if(i == 2) {
						dataOutput.writeInt(gTally[j]);

					} else if(i == 3) {
						dataOutput.writeInt(tTally[j]);

					}
				}
			}
			
			fileOutput.close();
			dataOutput.close();
		} catch (IOException e) {
			e.printStackTrace();
		}
    }

}