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
import dataStructure.ProtonFileHeader;

import net.sf.samtools.SAMFileHeader;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class CheKProton {

	/**
	 * @param args
	 * @throws IOException
	 */
	static long readCount = 0;

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
		TreeSet<IntervalAbs> intervalTSet = new TreeSet<IntervalAbs>();
		try (BufferedReader br = new BufferedReader(new FileReader(bedfile))) {
			// close by autoclosable
			String s;
			while ((s = br.readLine()) != null) {
				IntervalAbs intervalAbs = new IntervalAbs(s).setInterval();
				intervalTSet.add(intervalAbs);
			}

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
		endTime = Math.abs(System.nanoTime()) - startTime;
		System.err.println("IntervalList made ...." + intervalTSet.size()
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
		SAMFileHeader sfh =sfr.getFileHeader();
		ProtonFileHeader pfh = new ProtonFileHeader(sfh);
		for (SAMRecord samRecord : sfr) {
			readCount++;
			System.out.println(samRecord.getReadPairedFlag());
			String qname = samRecord.getReadName();
			int size = samRecord.getReadLength();
			int flags = samRecord.getFlags();
			String ref = samRecord.getReferenceName();
			String cig = samRecord.getCigarString();
			byte[] qual = samRecord.getBaseQualities();
			String qualString = samRecord.getBaseQualityString();
			int mQuality = samRecord.getMappingQuality();
			// ReferenceIndex seems to be the Index which is 0-based Chromosome
			// Number
			int readIndex = samRecord.getReferenceIndex();
			boolean isReadNeg = samRecord.getReadNegativeStrandFlag();
			boolean isReadUnmapped = samRecord.getReadUnmappedFlag();
			int alStart = samRecord.getAlignmentStart();
			int alEnd = samRecord.getAlignmentEnd();

			// get the offset for the absolute localization
			long offset = ChromosomeOffset.getChromosomeOffsetbyNumber(
					(short) (readIndex + 1)).getOffset();
			long absAlStart = 0;
			long absAlEnd = 0;

			if (!isReadNeg) {
				absAlStart = alStart + offset;
				absAlEnd = alEnd + offset;
			}
			if (isReadNeg) { // die Startkoordinate von reversen Reads ist im
								// BAM File nur über End-Koordinate minus Länge
								// des Reads aus dem Feld size + 1 zu berechnen
				absAlStart = offset + alEnd - size + 1;
				absAlEnd = offset + alEnd;
			}

			// sieht so aus, als müßte ich von dem Read ein IntervalAbs
			// erzeugen,
			// dass ich floor() übergeben kann, um das richtige
			// IntervalAbs- Objekt zurückzubekommen.

			IntervalAbs tempRead = new IntervalAbs((short) (readIndex + 1),
					absAlEnd, absAlStart, (int) (absAlEnd - absAlStart), null);

			IntervalAbs floorInterval = intervalTSet.floor(tempRead);
			endTime = Math.abs(System.nanoTime()) - startTime;
			if (floorInterval != null) {
				// System.out.println("Intervall: " +
				// floorInterval.getStartAbs()
				// + " - " + floorInterval.getEndAbs());
					System.out
							.printf("%.2f sec. %4d %-42s %5x %5s %8d - %-8d %4d %-14s %-8d bp (%d bp) %b%n",
									(double) (endTime / 1E9), readCount, qname,
									flags, ref, absAlStart, absAlEnd, mQuality,
									cig, (alEnd - alStart + 1), size, 
									isReadNeg);

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
			for (Iterator<IntervalAbs> iter = intervalTSet.iterator(); iter
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
			for (IntervalAbs interval : intervalTSet) {
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
				.printf("ReadCount: %,8d%n",
						readCount);
	} // end main
} // end class
