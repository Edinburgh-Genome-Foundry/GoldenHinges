from goldenhinges import OverhangsSelector
selector = OverhangsSelector(gc_min=0.25, gc_max=0.75,
	                         differences=2, time_limit=2)
collection = selector.generate_overhangs_set(n_overhangs=18)
print collection