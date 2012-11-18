package algorithm;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.Locale;
import java.util.TreeSet;

import dataStructure.ChromosomeOffset;
import dataStructure.FragmentReadEntry;
import dataStructure.PairedReadEntry;
import dataStructure.ReadEntry;

import net.sf.samtools.BAMRecord;
import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class CheKov {

	/**
	 * @param args
	 * @throws IOException
	 */
	private ReadEntry re;

	public static void main(String[] args) {
		long startTime = Math.abs(System.nanoTime());
		long endTime;
		String bedfile = args[0];
		String bamfile = args[1];
		String outfile = args[2];
		String missedBEDFile = args[3];
		// we need to do IMPORTANTLY some useful things with this parameter
		// Interval.INTERVAL_THRESHOLD = Integer.parseInt(args[4]);

		Locale.setDefault(Locale.US);
		TreeSet<IntervalAbs> intervalTreeSet = new TreeSet<IntervalAbs>();
		try (BufferedReader br = new BufferedReader(new FileReader(bedfile))) {
			// close by autoclosable
			String s;
			while ((s = br.readLine()) != null) {
				IntervalAbs intervalAbs = new IntervalAbs(s).setInterval();
				intervalTreeSet.add(intervalAbs);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		endTime = Math.abs(System.nanoTime()) - startTime;
		System.err.println("IntervalList made ...." + intervalTreeSet.size()
				+ " in " + (double) endTime / 1000000000 + " s");
		// TreeSet intervalTSet was filled with IntervalAbs Objects each
		// initialized with (short chr, long startAbs, long endAbs, int size)
		// AND ArrayList<Byte> coverage = new ArrayList<Byte>(size) each field
		// initialized with (byte)0;
		//

		/*
		 * Java 7 - Mehr als eine Insel von Christian Ullenboom Das Handbuch zu
		 * den Java SE-Bibliotheken: Die Klasse java.util.TreeSet implementiert
		 * ebenfalls wie HashSet die Set-Schnittstelle, verfolgt aber eine
		 * andere Implementierungsstrategie. Ein TreeSet verwaltet die Elemente
		 * immer sortiert (intern werden die Elemente in einem balancierten
		 * Binärbaum gehalten). Zugiffsmethoden sind ceiling(E e) = least
		 * element greater than or equal, higher(E e) = least element strictly
		 * greater than element E, floor(E e), lower(E e)
		 */

		/*
		 * d.h. eine eigene Implementierung eines BALANCIERTEN Binärbaums ist im
		 * Prinzip nicht nötig. Evtl. später Performance
		 */
		startTime = Math.abs(System.nanoTime());
		SAMFileReader sfr = new SAMFileReader(new File(bamfile));
		sfr.setValidationStringency(ValidationStringency.LENIENT);
		for (SAMRecord samRecord : sfr) {
			ReadEntry.readCount++;
			ReadEntry re = null;
			if (!samRecord.getReadPairedFlag()) {
				re = new FragmentReadEntry(samRecord);
			} else {
				re = new PairedReadEntry(samRecord);
			}
			re.analyzeCoverage();

			IntervalAbs floorInterval = intervalTreeSet.floor(tempRead);
			endTime = Math.abs(System.nanoTime()) - startTime;
			if (floorInterval != null) {

				if (floorInterval.readHitsAbsInterval(tempRead.getEndAbs(),
						tempRead.getStartAbs())) {
					// System.out.println(floorInterval.getCoverage());
					continue;
				}
			}

		} // end for
		/*
		 * here calculate missed areas in the intervals coverages. Strategy: 1.
		 * first iterate through the TreeSet 2. take the coverage ArrayList and
		 * iterate through it. 3. create a StringBuffer Instance, start
		 * appending when coverage at the position is 0 AND the reference to the
		 * StringBuffer is null. 4. Fill the second half of the StringBuffer if
		 * the bp-coverage != 0 again and the StringBuffer is != null. The
		 * StringBuffer is written to the FileWriter via a BufferdWriter, and
		 * the reference of the StringBuffer is nulled. This ensures that as
		 * many missed intervals can be identified in any BED interval.
		 */
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(
				missedBEDFile))) {
			for (Iterator<IntervalAbs> iter = intervalTreeSet.iterator(); iter
					.hasNext();) {
				IntervalAbs interval = iter.next();
				String chr = ChromosomeOffset.getChromosomeOffsetbyNumber(
						interval.getChr()).getChromosomeName();
				long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
						interval.getChr()).getOffset();
				int intervalStart = (int) (interval.getStartAbs() - offset);
				String geneName = interval.getGeneName();
				ArrayList<Integer> coverage = interval.getCoverage();
				StringBuffer zeroIntervall = null;

				for (int i = 0; i < coverage.size(); i++) {
					if (coverage.get(i) == 0 && zeroIntervall == null) {
						zeroIntervall = new StringBuffer(chr + "\t"
								+ (intervalStart + i) + "\t");
					}
					if (coverage.get(i) != 0 && zeroIntervall != null) {
						zeroIntervall.append((intervalStart + i) + "\t"
								+ geneName + "\n");
						bw.write(zeroIntervall.toString());
						zeroIntervall = null;
					}

				}

			}
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

		// Ausgabe Intervalle und coverages in Datei aus args[2] = outfile
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(outfile))) {
			for (IntervalAbs interval : intervalTreeSet) {
				long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
						interval.getChr()).getOffset();
				String line = String.format(
						"%s:%,d-%,d\t%s\t(%5d) %s%n",
						ChromosomeOffset.getChromosomeOffsetbyNumber(
								interval.getChr()).getChromosomeName(),
						(interval.getStartAbs() - offset),
						(interval.getEndAbs() - offset),
						interval.getGeneName(), interval.getSize(),
						interval.getCoverage());
				bw.write(line);
			}
		} catch (IOException ex) {
			ex.printStackTrace();
		}

		System.out
				.printf("ReadCount: %,8d%nReadMateUnmapped: %,8d%nReadMateChromosomeSplitted: %,8d%nBothSameOrientated: %,8d%nDist >1000: %,8d%n",
						readCount, readMateUnmappedCount,
						readMateSplittedChromosome, bothSameOrientated,
						highDistance);
	} // end main
} // end class
