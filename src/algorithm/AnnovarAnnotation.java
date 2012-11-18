package algorithm;

public class AnnovarAnnotation implements Comparable<Object> {
	private String annotation;
	private Integer count;
	

	public AnnovarAnnotation(String annotation) {
		this.setAnnotation(annotation);
		this.count = 1;
	}

	public String getAnnotation() {
		return annotation;
	}

	public void setAnnotation(String annotation) {
		this.annotation = annotation;

	}

	public Integer getCount() {
		return count;
	}

	public void setCount(Integer count) {
		this.count = count;
	}

	@Override
	public boolean equals(Object obj) {
		if (obj == null)
			return false;
		if (obj == this)
			return true;

		if (!(obj instanceof AnnovarAnnotation))
			return false;
		AnnovarAnnotation that = (AnnovarAnnotation) obj;
		return (that.annotation.equals(this.annotation));
	}

	@Override
	public int compareTo(Object obj) {
		AnnovarAnnotation that = (AnnovarAnnotation) obj;
		return this.annotation.compareTo(that.annotation);
	}
}
