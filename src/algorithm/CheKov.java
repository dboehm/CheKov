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

import logger.LogReadConfigExample;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;
import org.apache.commons.cli.PosixParser;
import org.apache.log4j.Logger;
import org.apache.log4j.PropertyConfigurator;

import dataStructure.ChromosomeOffset;
import dataStructure.FragmentReadEntry;
import dataStructure.IntervalAbs;
import dataStructure.PairedReadEntry;
import dataStructure.ReadEntry;
import dataStructure.TargetNucleotidePositionEntry;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;

public class CheKov {

	private static final long FIVE_SECS = 5000;
	/**
	 * @param args
	 * @throws IOException
	 */
	private static TreeSet<IntervalAbs> intervalTreeSet;
	private static long allBases = 0;
	private static long allHardClippedBases = 0;
	private static long allSoftClippedBases = 0;
	private static long allDeletedTaggedBases = 0;
	private static long allInsertedTaggedBases = 0;
	private static long allDeletedTaggedReads = 0;
	private static long allInsertedTaggedReads = 0;
	private static long allEitherDeletedOrInsertedTaggedReads = 0;
	private static long medianEffReadLength = 0;
	private static long medianRawReadLength = 0;
	private static long avEffReadLength = 0;
	private static long avRawReadLength = 0;
	private static int allReadsWithoutQuality = 0;
	private static TreeSet<TargetNucleotidePositionEntry> alteredNucleotidePositionsEntries = new TreeSet<>();
	private static long startTime = 0;
	private static long endTime = 0;
	private static IndexedFastaSequenceFile indexedFastaSequenceFile_Ref = null;

	public static void main(String[] args) {
		final Logger logger = Logger.getLogger(CheKov.class);
		PropertyConfigurator.configureAndWatch("config/log4j.properties",
				FIVE_SECS);
		logger.info("CheKovLogger started");
		String bedfile = null;
		String bamfile = null;
		String outfile = null;
		String missedBEDFile = null;
		// we need to do IMPORTANTLY some useful things with this parameter
		IntervalAbs.INTERVAL_THRESHOLD = 0;
		String refFile = null;
		String homopolymerFile = null;
		final Options allowedOptions = new Options();
		// Add all options from enum AppParameterCLI
		for (AppParameterCLI parameter : AppParameterCLI.values()) {
			allowedOptions.addOption(parameter.shortName, parameter.name,
					parameter.hasArgs, parameter.helpText);
		}
		final CommandLineParser parser = new PosixParser();
		try {
			final CommandLine commandLine = parser.parse(allowedOptions, args);
			homopolymerFile = commandLine.getOptionValue('h');
			refFile = commandLine.getOptionValue('r');
			IntervalAbs.INTERVAL_THRESHOLD = Integer.parseInt(commandLine
					.getOptionValue('t'));
			missedBEDFile = commandLine.getOptionValue('m');
			outfile = commandLine.getOptionValue('o');
			bedfile = commandLine.getOptionValue('b');
			bamfile = commandLine.getOptionValue('a');

		} catch (ParseException e3) {
			// TODO Auto-generated catch block
			HelpFormatter help = new HelpFormatter();
			help.printHelp("CheKov", allowedOptions);
			System.exit(0);
		}
		logger.info("BAM-File: '" + bamfile + "'");
		logger.info("BED-File: '" + bedfile + "'");
		logger.info("Reference-File: '" + refFile + "'");
		logger.info("Outfile: '" + outfile + "'");
		logger.info("missed-File: '" + missedBEDFile + "'");
		logger.info("Homopolymer-File: '" + homopolymerFile + "'");
		CheKov.setStartTime(Math.abs(System.nanoTime()));

		/*
		 * TreeSet intervalTreeSet is filled with IntervalAbs Objects each
		 * initialized with (short chr, long startAbs, long endAbs, int size)
		 * AND ArrayList<Integer> coverage = new ArrayList<>(size) each field
		 * initialized with (Integer) 0;
		 */
		Locale.setDefault(Locale.US);
		intervalTreeSet = new TreeSet<IntervalAbs>();
		try (BufferedReader br = new BufferedReader(new FileReader(bedfile))) {
			// closed by autoclosable
			String s;
			while ((s = br.readLine()) != null) {
				IntervalAbs intervalAbs = new IntervalAbs(s).setInterval();
				intervalTreeSet.add(intervalAbs);
			}

		} catch (FileNotFoundException ex) {
			ex.printStackTrace();
		} catch (IOException ex) {
			ex.printStackTrace();
		}
		CheKov.setEndTime(Math.abs(System.nanoTime()) - CheKov.getStartTime());
		logger.info("IntervalList made ...." + intervalTreeSet.size() + " in "
				+ (double) endTime / 1000000000 + " s");

		// Now start the ReferenceReader

		try {
			indexedFastaSequenceFile_Ref = new IndexedFastaSequenceFile(
					new File(refFile));
			if (!indexedFastaSequenceFile_Ref.isIndexed())
				System.exit(0);
		} catch (FileNotFoundException e2) {
			// TODO Auto-generated catch block
			e2.printStackTrace();
		}
		/*
		 * Java 7 - Mehr als eine Insel von Christian Ullenboom Das Handbuch zu
		 * den Java SE-Bibliotheken: " ... Ein TreeSet verwaltet die Elemente
		 * immer sortiert (intern werden die Elemente in einem balancierten
		 * Binärbaum gehalten). Zugiffsmethoden sind ceiling(E e) = least
		 * element greater than or equal, higher(E e) = least element strictly
		 * greater than element E, floor(E e), lower(E e) d.h. eine eigene
		 * Implementierung eines BALANCIERTEN Binärbaums ist im Prinzip nicht
		 * nötig. Evtl. später Performance
		 */
		CheKov.setStartTime(Math.abs(System.nanoTime()));
		SAMFileReader samFileReader = new SAMFileReader(new File(bamfile));
		samFileReader.setValidationStringency(ValidationStringency.LENIENT);
		for (SAMRecord samRecord : samFileReader) {
			// this is the count for all reads coming in
			if (samRecord.getReferenceName() == "chrM")
				continue;

			// check if length of read and length of quality array is same
			if (samRecord.getReadBases().length != samRecord.getBaseQualities().length) {
				CheKov.incrementReadsWithoutQualities();
				continue;
			}
			/*
			 * the result is: mapping coordinates getAlignmentStart() and
			 * getAlignmentEnd() have 1. have included the deleted Positions
			 * from the read as a gap. means from Cigar, if 1D is in Cigar, the
			 * length of Alignment is each 1 position added. (==> if you give
			 * the Byte[] of the Read and/or Quality AND the Cigar, then 1D
			 * leaves a gap in the Interval, and does not set coverage and
			 * quality for that position in the IntervalAbs. BUT coverage and
			 * quality need to be set in the TargetNucleotidePositionEntry) 2.
			 * have NOT included the soft-clipped Positions. In this analysis
			 * they do not contribute to coverage and quality. Coming from
			 * Byte[] of the Read and of quality these positions need be skipped
			 * 3. Insertions in Cigar (e.g. 1I) do contribute to coverage and
			 * quality of the TargetNucleotidePositionEntry, but NOT on the
			 * IntervalAbs.
			 * 
			 * If things are unclear uncomment the below output and check
			 */
			// check if length and mapping coordinates fit together
			// printSamRecord(samRecord);

			// check the output in detail
			// for (int i = 0; i < samRecord.getBaseQualities().length; i++)
			// System.out.print((char) samRecord.getReadBases()[i] + ""
			// + samRecord.getBaseQualities()[i] + " ");
			// System.out.println();
			ReadEntry.setReadCount(ReadEntry.getReadCount() + 1);
			ReadEntry readEntry = null;
			// if the Read is initially a Single Fragment Read, initialize a
			// FragmentReadEntry and count the numbers
			if (!samRecord.getReadPairedFlag()) {
				FragmentReadEntry.setFragmentReadCount(FragmentReadEntry
						.getFragmentReadCount() + 1);
				readEntry = new FragmentReadEntry(samRecord);
				// else if the Read is initially a Paired Read, initialize a
				// PairedReadEntry and count the numbers
			} else {
				PairedReadEntry.setPairedEndReadCount(PairedReadEntry
						.getPairedEndReadCount() + 1);
				readEntry = new PairedReadEntry(samRecord);
			}
			// this is a simple form of a progress bar, printing a line every 1
			// Million Reads
			if (ReadEntry.getReadCount() % 1_000_00 == 0) {
				CheKov.setEndTime(Math.abs(System.nanoTime())
						- CheKov.getStartTime());
				logger.info(String.format("%,d %,d %,d %5.2f %s  %,d %s",
						ReadEntry.getReadCount(),
						FragmentReadEntry.getFragmentReadCount(),
						PairedReadEntry.getPairedEndReadCount(),
						(double) CheKov.getEndTime() / 1000000000, "s",
						alteredNucleotidePositionsEntries.size(),
						"alteredNucleotidePositionsEntries"));
				CheKov.setStartTime(Math.abs(System.nanoTime()));
			}
			// calculate the coverage of each read on the targets represented by
			// intervalTreeSet
			readEntry.analyzeCoverage();
			// check the quality of the reads
			// readEntry.analyseQuality();
		} // end for
		CheKov.setEndTime(Math.abs(System.nanoTime()) - CheKov.getStartTime());

		/*
		 * get the two ReadLengthCounterArray for median calculation use static
		 * method medianFromArray to find the middle value
		 */
		CheKov.setMedianRawReadLength(CheKov.medianFromArray(FragmentReadEntry
				.getRawLengthCounterArray()));
		CheKov.setMedianEffReadLength(CheKov.medianFromArray(FragmentReadEntry
				.getEffLengthCounterArray()));
		/*
		 * get the coverage for each TargetNucleotidePositionEntry from TreeSet
		 * alteredNucleotidePositionsEntries just to show fast the result as
		 * output! need to be implemented in other way soon
		 */
		try (BufferedWriter bw = new BufferedWriter(new FileWriter(
				homopolymerFile))) {
			for (TargetNucleotidePositionEntry tnpe : alteredNucleotidePositionsEntries) {
				long positionInInterval = -1;
				IntervalAbs tempRead = new IntervalAbs((short) 0,
						tnpe.getPosInAbsGenome(), tnpe.getPosInAbsGenome(), 1,
						null);
				IntervalAbs floorInterval = CheKov.getIntervalTreeSet().floor(
						tempRead);
				if (floorInterval != null) {

					boolean hitted = false;
					ArrayList<Integer> cov = floorInterval.getCoverage();
					for (long i = floorInterval.getStartAbs(); i < floorInterval
							.getEndAbs(); i++) {
						if (i >= tempRead.getEndAbs()
								&& i <= tempRead.getStartAbs()) {
							positionInInterval = tnpe.getPosInAbsGenome()
									- floorInterval.getStartAbs();
							tnpe.setCoverage(cov.get((int) positionInInterval));
							hitted = true;
						} // end inner if
					} // end inner for

					if (hitted) {
						bw.write(String.format("%s%n", tnpe));
					}

				} // end outer if
			} // end outer for
		} catch (IOException e1) {
			// TODO Auto-generated catch block
			e1.printStackTrace();
		}

		// finale Ausgabe
		System.out.println(alteredNucleotidePositionsEntries.size());

		printQCResult();

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

	}// end main

	public static void printSamRecord(SAMRecord samRecord) {
		System.out
				.println("\n"
						+ samRecord.getReadName()
						+ " "
						+ samRecord.getAlignmentStart()
						+ "-"
						+ samRecord.getAlignmentEnd()
						+ "  AlignmentLength: "
						+ (samRecord.getAlignmentEnd()
								- samRecord.getAlignmentStart() + 1)
						+ " ReadLength: " + (samRecord.getReadBases().length)
						+ " : " + samRecord.getCigarString());
	}

	public static void printQCResult() {
		System.out.printf("%-35s%,12d%n%-35s%,12d%n%-35s%,12d%n",
				"All Reads: ", ReadEntry.getReadCount(), "FragmentReads:",
				FragmentReadEntry.getFragmentReadCount(), "PairedReads:",
				PairedReadEntry.getPairedEndReadCount());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n", "OnTargetReadCounts: ",
				FragmentReadEntry.getOnTargetReadCount(),
				"OffTargetReadCount:",
				FragmentReadEntry.getOffTargetReadCount());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n%-35s%,12d%n", "Seq-bp",
				CheKov.getAllBases(), "SoftCipped-bp:",
				CheKov.getAllSoftClippedBases(), "HardClipped-bp:",
				CheKov.getAllHardClippedBases());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n%-35s%,12d%n",
				"Unmapped Reads", FragmentReadEntry.getReadUnmappedCount(),
				"Unmapped Mate:", PairedReadEntry.getReadMateUnmappedCount(),
				"Mate Splitted Chr:",
				PairedReadEntry.getReadMateSplittedChromosome());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n", "Pairs same orientation",
				PairedReadEntry.getBothSameOrientated(),
				"Pair distance > 1000 bp:", PairedReadEntry.getHighDistance());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n", "Deleted Tagged bp:",
				CheKov.getAllDeletedTaggedBases(), "Inserted Tagged bp:",
				CheKov.getAllInsertedTaggedBases());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n%-35s%,12d%n",
				"Deleted Tagged Reads:", CheKov.getAllDeletedTaggedReads(),
				"Inserted Tagged Reads:", CheKov.getAllInsertedTaggedReads(),
				"Deleted- OR Inserted Tagged Reads:",
				CheKov.getAllEitherDeletedOrInsertedTaggedReads());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n", "AvRawReadLength:",
				CheKov.getAvRawReadLength() / ReadEntry.getReadCount(),
				"AvEffReadLength:",
				CheKov.getAvEffReadLength() / ReadEntry.getReadCount());
		System.out.printf("%-35s%,12d%n%-35s%,12d%n", "MedianRawReadLength:",
				CheKov.getMedianRawReadLength(), "MedianEffReadLength:",
				CheKov.getMedianEffReadLength());
	}

	public static int medianFromArray(int[] array) {
		int counter = ReadEntry.getReadCount();
		for (int i = 0; i < 2000; i++) {
			counter = counter - array[i];
			if (counter <= (int) ReadEntry.getReadCount() / 2) {
				return i;
			}
		}
		return 0;
	}

	// Getter and Setter
	public static TreeSet<IntervalAbs> getIntervalTreeSet() {
		return intervalTreeSet;
	}

	public static void setIntervalTreeSet(TreeSet<IntervalAbs> intervalTreeSet) {
		CheKov.intervalTreeSet = intervalTreeSet;
	}

	public static long getAllBases() {
		return allBases;
	}

	public static void setAllBases(long allBases) {
		CheKov.allBases = allBases;
	}

	public static long getAllHardClippedBases() {
		return allHardClippedBases;
	}

	public static void setAllHardClippedBases(long allHardClippedBases) {
		CheKov.allHardClippedBases = allHardClippedBases;
	}

	public static long getAllSoftClippedBases() {
		return allSoftClippedBases;
	}

	public static void setAllSoftClippedBases(long allSoftClippedBases) {
		CheKov.allSoftClippedBases = allSoftClippedBases;
	}

	public static long getAllDeletedTaggedBases() {
		return allDeletedTaggedBases;
	}

	public static void setAllDeletedTaggedBases(long allDeletedTaggedBases) {
		CheKov.allDeletedTaggedBases = allDeletedTaggedBases;
	}

	public static long getAllInsertedTaggedBases() {
		return allInsertedTaggedBases;
	}

	public static void setAllInsertedTaggedBases(long allInsertedTaggedBases) {
		CheKov.allInsertedTaggedBases = allInsertedTaggedBases;
	}

	public static long getAllDeletedTaggedReads() {
		return allDeletedTaggedReads;
	}

	public static void setAllDeletedTaggedReads(long allDeletedTaggedReads) {
		CheKov.allDeletedTaggedReads = allDeletedTaggedReads;
	}

	public static long getAllInsertedTaggedReads() {
		return allInsertedTaggedReads;
	}

	public static void setAllInsertedTaggedReads(long allInsertedTaggedReads) {
		CheKov.allInsertedTaggedReads = allInsertedTaggedReads;
	}

	public static long getAllEitherDeletedOrInsertedTaggedReads() {
		return allEitherDeletedOrInsertedTaggedReads;
	}

	public static void setAllEitherDeletedOrInsertedTaggedReads(
			long allEitherDeletedOrInsertedTaggedReads) {
		CheKov.allEitherDeletedOrInsertedTaggedReads = allEitherDeletedOrInsertedTaggedReads;
	}

	public static long getMedianEffReadLength() {
		return medianEffReadLength;
	}

	public static void setMedianEffReadLength(long medianEffReadLength) {
		CheKov.medianEffReadLength = medianEffReadLength;
	}

	public static long getMedianRawReadLength() {
		return medianRawReadLength;
	}

	public static void setMedianRawReadLength(long medianRawReadLength) {
		CheKov.medianRawReadLength = medianRawReadLength;
	}

	public static long getAvEffReadLength() {
		return avEffReadLength;
	}

	public static void setAvEffReadLength(long avEffReadLength) {
		CheKov.avEffReadLength = avEffReadLength;
	}

	public static long getAvRawReadLength() {
		return avRawReadLength;
	}

	public static void setAvRawReadLength(long avRawReadLength) {
		CheKov.avRawReadLength = avRawReadLength;
	}

	public static TreeSet<TargetNucleotidePositionEntry> getAlteredNucleotidePositionsEntries() {
		return alteredNucleotidePositionsEntries;
	}

	public static void setAlteredNucleotidePositionsEntries(
			TreeSet<TargetNucleotidePositionEntry> alteredNucleotidePositionsEntries) {
		CheKov.alteredNucleotidePositionsEntries = alteredNucleotidePositionsEntries;
	}

	public static long getStartTime() {
		return startTime;
	}

	public static void setStartTime(long startTime) {
		CheKov.startTime = startTime;
	}

	public static long getEndTime() {
		return endTime;
	}

	public static void setEndTime(long endTime) {
		CheKov.endTime = endTime;
	}

	public static IndexedFastaSequenceFile getIndexedFastaSequenceFile_Ref() {
		return indexedFastaSequenceFile_Ref;
	}

	public static void setIndexedFastaSequenceFile_Ref(
			IndexedFastaSequenceFile indexedFastaSequenceFile_Ref) {
		CheKov.indexedFastaSequenceFile_Ref = indexedFastaSequenceFile_Ref;
	}

	public static int getAllReadsWithoutQuality() {
		return allReadsWithoutQuality;
	}

	public static void incrementReadsWithoutQualities() {
		CheKov.allReadsWithoutQuality++;
	}

} // end class
