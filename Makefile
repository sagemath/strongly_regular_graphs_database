database:
	./website_to_tmp.sh
	python tmp_to_db.py
	rm -f brouwer.tmp
clean:
	rm -f brouwer_srg_database.txt
	rm -f brouwer_srg_database.json
	rm -f brouwer.tmp

