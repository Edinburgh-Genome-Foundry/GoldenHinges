from goldenhinges import OverhangsSelector
selector = OverhangsSelector(gc_min=0.25, gc_max=0.5,
	                         differences=2, time_limit=1)
collection = selector.generate_overhangs_set(n_cliques=100)
collection = selector.generate_overhangs_set(start_at=len(collection))
print ("Found %d overhangs: " %len(collection), collection)
