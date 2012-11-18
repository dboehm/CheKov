package algorithm;

import java.util.ArrayList;

import dataStructure.ChromosomeOffset;

public class IntervalAbs implements Comparable<IntervalAbs> {

	private String line;
	private short chr;
	private long intervalStartAbs;
	private long intervalEndAbs;
	private String geneName;
	private int intervalSize;
	private ArrayList<Integer> coverage;

	public IntervalAbs(String line) {
		this.line = line;
	}

	public IntervalAbs(short chr, long startAbs, long endAbs, int size, String geneName) {
		this.chr = chr;
		this.intervalStartAbs = startAbs;
		this.intervalEndAbs = endAbs;
		this.intervalSize = size;
		this.geneName = geneName;
		this.coverage = new ArrayList<Integer>(size);
		for (int i = 0; i < size; i++) {
			this.coverage.add(0);
		}
	}

	public IntervalAbs setInterval() {
		String[] lines = this.line.split("\t");
		short chr = ChromosomeOffset.chromosomeNumber(lines[0]);
		long startAbs = Long.parseLong(lines[1])
				+ ChromosomeOffset.getChromosomeOffsetbyNumber(chr).getOffset();
		long endAbs = Long.parseLong(lines[2])
				+ ChromosomeOffset.getChromosomeOffsetbyNumber(chr).getOffset();
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

	
	public String getGeneName() {
		return geneName;
	}

	public void setGeneName(String geneName) {
		this.geneName = geneName;
	}

	@Override
	public int compareTo(IntervalAbs o) {
		if (this.equals(o))
			return 0;
		// the range is from ~3E9 which is to much for an int, /10 maximizes
		// values
		return (int) (this.intervalStartAbs / 10 - o.intervalStartAbs / 10);
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

	public boolean readHitsAbsInterval(long readAbsStart, long readAbsEnd) {
		boolean hitted = false;
		ArrayList<Integer> cov = this.getCoverage();
		for (long i = this.intervalStartAbs; i < this.intervalEndAbs; i++) {
			if (i >= readAbsStart && i < readAbsEnd) {
				cov.set((int) (i - this.intervalStartAbs),
						cov.get((int) (i - this.intervalStartAbs)) + 1);
				hitted = true;
			}
		}
		this.setCoverage(cov);
		return hitted;
	}
}