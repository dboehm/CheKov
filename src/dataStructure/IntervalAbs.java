package dataStructure;

import java.util.ArrayList;

import algorithm.CheKov;

public class IntervalAbs implements Comparable<IntervalAbs> {

	public static int INTERVAL_THRESHOLD;
	private String line;
	private short chr;
	private long intervalStartAbs;
	private long intervalEndAbs;
	private String geneName;
	private int intervalSize;
	private ArrayList<Integer> coverage;
	private ArrayList<Integer> quality;

	public IntervalAbs(String line) {
		this.line = line;
	}

	public IntervalAbs(short chr, long startAbs, long endAbs, int size,
			String geneName) {
		this.chr = chr;
		this.intervalStartAbs = startAbs;
		this.intervalEndAbs = endAbs;
		this.intervalSize = size;
		this.geneName = geneName;
		// fill ArrayList coverage with 0
		this.coverage = new ArrayList<Integer>(size);
		for (int i = 0; i < size; i++) {
			this.coverage.add(0);
		}
		// fill ArrayList quality with 0
		this.quality = new ArrayList<Integer>(size);
		for (int i = 0; i < size; i++) {
			this.quality.add(0);
		}
	}

	public IntervalAbs setInterval() {
		String[] lines = this.line.split("\t");
		short chr = ChromosomeOffset.chromosomeNumber(lines[0]);
		long startAbs = Long.parseLong(lines[1])
				+ ChromosomeOffset.getChromosomeOffsetbyNumber(chr).getOffset()
				- INTERVAL_THRESHOLD;
		long endAbs = Long.parseLong(lines[2])
				+ ChromosomeOffset.getChromosomeOffsetbyNumber(chr).getOffset()
				+ INTERVAL_THRESHOLD;
		intervalSize = (int) (endAbs - startAbs);
		String[] fields = lines[3].split(":");

		return new IntervalAbs(chr, startAbs, endAbs, intervalSize, fields[2]);
	}

	public String getLine() {
		return line;
	}

	public void setLine(String line) {
		this.line = line;
	}

	public long getStartAbs() {
		return intervalStartAbs;
	}

	public void setStartAbs(long startAbs) {
		this.intervalStartAbs = startAbs;
	}

	public long getEndAbs() {
		return intervalEndAbs;
	}

	public void setEndAbs(long endAbs) {
		this.intervalEndAbs = endAbs;
	}

	public int getSize() {
		return intervalSize;
	}

	public void setSize(int size) {
		this.intervalSize = size;
	}

	public short getChr() {
		return chr;
	}

	public void setChr(short chr) {
		this.chr = chr;
	}

	public ArrayList<Integer> getCoverage() {
		return coverage;
	}

	public void setCoverage(ArrayList<Integer> coverage) {
		this.coverage = coverage;
	}

	public ArrayList<Integer> getQuality() {
		return quality;
	}

	public void setQuality(ArrayList<Integer> quality) {
		this.quality = quality;
	}

	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	@Override
	public int compareTo(IntervalAbs other) {
		if (this.equals(other))
			return 0;
		// the range is from ~3E9 which is to much for an int, /10 maximizes
		// values
		return (int) (this.intervalStartAbs / 2 - other.intervalStartAbs / 2);
		// the before implementation seems to be buggy, this following
		// implementation even more :-(
		// return (this.getChr()
		// * ((int) (this.intervalStartAbs - ChromosomeOffset
		// .getChromosomeOffsetbyNumber(this.getChr()).getOffset())))
		// - ((other.getChr()
		// * ((int) (other.intervalStartAbs - ChromosomeOffset
		// .getChromosomeOffsetbyNumber(other.getChr()).getOffset()))));
	}

	@Override
	public boolean equals(Object obj) {
		// equals in terms of absolute intervals means same start and end,
		// ToDo: What is about overlaps ??
		if (obj instanceof IntervalAbs) {
			IntervalAbs new_name = (IntervalAbs) obj;
			if (this.intervalStartAbs == new_name.intervalStartAbs)
				return true;
		}
		return super.equals(obj);
	}

	public boolean readHitsAbsInterval(long readAbsStart, long readAbsEnd,
			byte[] qualities, ArrayList<String> cigarTokens) {
		boolean hitted = false;
		// get the copy of the so far ArrayLists of the IntervalAbs
		ArrayList<Integer> cov = this.getCoverage();
		ArrayList<Integer> qual = this.getQuality();
		int relInReadPosition = 0;
		int relInIntervalPosition = 0;
		int inCigarPosition = 0;
		/*
		 * there are 3 + 1 scenarios: case 1: Read starts BEFORE the Interval
		 * and ends IN the Interval or following the Interval ==>
		 * intervalStartAbs >= readAbsStart && intervalStartAbs <= readAbsEnd
		 * case 2: Read starts ON the Interval (either equal or higher position)
		 * intervalStartAbs < readAbsStart && intervalEndAbs >= readAbsStart ==>
		 * relInReadPosition = 0; relInIntervalPosition = readAbsStart -
		 * this.intervalStartAbs case 3: ==> read does not overlap with Interval
		 * ==> hitted = false intervalStartAbs < readAbsStart && intervalEndAbs
		 * < readAbsStart case 4: seems to be a bug, perhaps because of the
		 * implementation of the compareTo() - method in IntervalAbs, we changed
		 * the implementation of the compareTo method to divide 2 and removed
		 * bug case 4, which was IntervalAbs was at all bigger than ReadInterval
		 * without overlap - still checked by assert
		 */
		if (this.intervalStartAbs > readAbsStart
				&& this.intervalStartAbs <= readAbsEnd + 1) {
			relInReadPosition = (int) (this.intervalStartAbs - readAbsStart + 1);
			relInIntervalPosition = 0;
			hitted = true;
		} else if (this.intervalStartAbs <= readAbsStart
				&& this.intervalEndAbs > readAbsStart) {
			relInReadPosition = 0;
			relInIntervalPosition = (int) (readAbsStart - this.intervalStartAbs);
			hitted = true;
		} else if (this.intervalStartAbs <= readAbsStart
				&& this.intervalEndAbs <= readAbsStart) {
			return hitted; // it is false here
		} else {
			assert !(this.intervalStartAbs > readAbsStart && this.intervalEndAbs > readAbsStart) : "Return interval of floor() completely downstream of read. Check compareTo()";
			return hitted; // it is false here
		}
		System.out.println("In Interval " + relInIntervalPosition
				+ "  InRead : " + relInReadPosition);
		for (int i = 0; i < cigarTokens.size(); i += 2) {
			int tokenCount = Integer.parseInt(cigarTokens.get(i));
			String tokenType = cigarTokens.get(i + 1);
			int oldRelInReadPosition = inCigarPosition;
			switch (tokenType) {
			case "S":
				inCigarPosition += tokenCount;
				break;

			case "M":
				inCigarPosition += tokenCount;
				break;
			case "D":
				for (int i1 = 1; i1 <= tokenCount; i1++)
					System.out.print("D");
				break;
			case "I":
				inCigarPosition += tokenCount;
				break;
			}
			for (int j = oldRelInReadPosition; j < inCigarPosition; j++)
				System.out.print(qualities[j] + "+");
			System.out.print("+++");
		}
		System.out.println();
		// this prints the complete quality byte[] array
		for (int j = 0; j < qualities.length; j++)
			System.out.print(qualities[j] + ":");
		System.out.println();
		System.out.println(" " + qualities.length + " readAbsStart-ReadAbsEnd="
				+ readAbsStart + "-" + readAbsEnd + "="
				+ (readAbsEnd - readAbsStart + 1));
		// this was the old implementation
		// for (long i = this.intervalStartAbs; i < this.intervalEndAbs; i++) {
		// if (i >= readAbsStart && i < readAbsEnd) {
		// relInIntervalPosition = (int) (i - this.intervalStartAbs);
		// relInReadPosition = (int) (i - readAbsStart);
		// cov.set((int) (i - this.intervalStartAbs),
		// cov.get((int) (i - this.intervalStartAbs)) + 1);
		// if (qualities.length == 0) {
		// CheKov.incrementReadsWithoutQualities();
		// } else {
		// qual.set((int) (i - this.intervalStartAbs),
		// qual.get((int) (i - this.intervalStartAbs)) + 1);
		// }
		// hitted = true;
		// }
		// }
		this.setCoverage(cov);
		return hitted;
	}

	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return String.format("%s%d:%,d-%,d (%,d bp)", "chr", this.getChr(),
				(this.getStartAbs() - ChromosomeOffset.offset(getChr())),
				(this.getEndAbs() - ChromosomeOffset.offset(getChr())),
				this.getSize());
	}
}