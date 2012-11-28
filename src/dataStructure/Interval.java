package dataStructure;

import java.util.ArrayList;

public class Interval {
	private int start, end;
	private String chr;
	private String line;
	private int size;
	private ArrayList<Byte> coverage;

	public Interval(String chr, int start, int end, int size) {
		this.chr = chr;
		this.start = start;
		this.end = end;
		this.size = size;
		System.out.println("FILL");
		this.coverage = new ArrayList<Byte>(size);
		for (int i = 0; i < size; i++) {
			System.out.println("FILL");
			this.coverage.add((byte) 0);
		}
	}

	public Interval(String line) {
		this.line = line;
	}

	public Interval setInterval() {
		String[] lines = this.line.split("\t");
		int start = Integer.parseInt(lines[1]);
		int end = Integer.parseInt(lines[2]);

		return new Interval(lines[0].toString(), start, end, (end - start));

	}

	public int getStart() {
		return start;
	}

	public void setStart(int start) {
		this.start = start;
	}

	public int getEnd() {
		return end;
	}

	public void setEnd(int end) {
		this.end = end;
	}

	public String getChr() {
		return chr;
	}

	public void setChr(String chr) {
		this.chr = chr;
	}

	public String getLine() {
		return line;
	}

	public void setLine(String line) {
		this.line = line;
	}

	public int getSize() {
		return size;
	}

	public void setSize(int size) {
		this.size = size;
	}

	public ArrayList<Byte> getCoverage() {
		return coverage;
	}

	public void setCoverage(ArrayList<Byte> coverage) {
		this.coverage = coverage;
	}

}
