package algorithm;

import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Set;
import java.util.TreeSet;

import dataStructure.ChromosomeOffset;
import dataStructure.GATK_SNV_entry;

public class ChekAnnovar {
	private static Set<AnnovarAnnotation> annotationEntries = new TreeSet<>();
	public static int count = 0;
	private static List<Annovar_GATK_SNV_Entry> annovarGATKSNV_EntryList = new ArrayList<>();

	/**
	 * @param args
	 */
	public static void main(String[] args) {

		try (BufferedReader br = new BufferedReader(new FileReader(args[0]))) {
			String line;
			while ((line = br.readLine()) != null) {
				String[] lineFields = line.split("\t");
				createAnnotationSet(lineFields[0]);
				Annovar_GATK_SNV_Entry annovarGATKSNV_Entry = new Annovar_GATK_SNV_Entry(
						lineFields[0], lineFields[1],
						new GATK_SNV_entry(lineFields[7], Long
								.parseLong(lineFields[8]), Long
								.parseLong(lineFields[8])
								+ ChromosomeOffset.offset(ChromosomeOffset
										.chromosomeNumber(lineFields[7])),
								ChromosomeOffset
										.chromosomeNumber(lineFields[7]),
								lineFields[9], lineFields[10], lineFields[11],
								Double.parseDouble(lineFields[12]),
								lineFields[13], lineFields[14], lineFields[15],
								lineFields[16]));
				annovarGATKSNV_EntryList.add(annovarGATKSNV_Entry);
			}


//			for (Annovar_GATK_SNV_Entry annoGATK_Entry : annovarGATKSNV_EntryList) {
//				System.out.printf("%-20s %-10s %-5s %5d %10d %12d %-12s%n",
//						annoGATK_Entry.getAnnotationType(), annoGATK_Entry
//								.getGene(), annoGATK_Entry.getChr(),
//						annoGATK_Entry.getChr_nr(), annoGATK_Entry
//								.getGenomPosition(), annoGATK_Entry
//								.getAbsGenomPosition(), annoGATK_Entry
//								.getGatk_SNV_entry().getDbSNP_entry());
//			}

			for (AnnovarAnnotation anno : annotationEntries) {
				System.out.printf("%-20s%10d%n", anno.getAnnotation(),
						anno.getCount());
				count += anno.getCount();
			}
			System.out.printf("%-20s%10d%n", "Total", count);

		} catch (FileNotFoundException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private static void createAnnotationSet(String s) {
		AnnovarAnnotation newAnno = new AnnovarAnnotation(s);
		if (!annotationEntries.contains(newAnno)) {
			annotationEntries.add(newAnno);
		} else {
			for (AnnovarAnnotation anno : annotationEntries) {
				if (anno.equals(newAnno)) {
					anno.setCount(anno.getCount() + 1);
				}
			}
		}
	}
} // end class
