// to change track name width -> in ../browser/sukn.js, browser.leftColumnWidth=210;
var _DEBUG_ = false;

var tmp_json_url 	= "https://egg2.wustl.edu/web_portal_cache/";
//var tmp_json_url 	= "https://mitra.stanford.edu/kundaje/leepc12/web_portal_cache/";

var host_url 		= "https://epigenomegateway.wustl.edu";
var wubrowser_url 	= "https://epigenomegateway.wustl.edu/browser/"

var track_num_alert	= 500; // warning for more than # tracks

var width_trackname 	= 150;

var height_row_bigwig_default	= 30;
var height_row_cat_default	= 10; 
var height_row_qcat_default 	= 150;
var height_row_catmat_default 	= 150;
var height_row_bed_default	= 30; 
var height_row_interact_default	= 120; 
var height_row_hammock_default	= 18;

// extra columns for grid visulation; chr state models 15(primary), 18(expanded) and 25(imputed)
var assay_num_vis = 3;
var assay_name_chr_prim15  	= "15-STATE MODEL";
var assay_name_chr_exp18   	= "18-STATE MODEL";
var assay_name_imp_chr_12m 	= "25-STATE MODEL (imputed)";
var assay_DNAMethyl 		= "DNAMethylSBS";

// GWAS analysis thread
var thread_sleep_ms = 5000; // 
var curr_job_id_str = ""; // current job id , if empty, there will be no heartbeat CGI to check GWAS log

function debug_log( str )
{
	if (_DEBUG_) console.log( str );
}

function component_to_hex( c ) 
{
	var hex = c.toString(16);
	return hex.length == 1 ? "0" + hex : hex;
}

function rgb_to_hex( rgb ) 
{
	var arr = rgb.split(",");	
	var r = parseInt( arr[0] );
	var g = parseInt( arr[1] );
	var b = parseInt( arr[2] );
	return "#" + component_to_hex(r) + component_to_hex(g) + component_to_hex(b);
}

function embed_wubrowser( id, json_url, genome )
{
	var elem = document.getElementById( id );
	elem.style.display = "inline-block";

	embed_washugb({
		leftSpaceWidth:200,
		host:host_url,
		container:elem,
		genome:genome,
		maxTrackHeight:100,
		datahub:json_url,
	});
}

function visualize_data( data, filter, embed_id, sort, genome ) // cgi-bin from visualize_data.c, filter is array
{
	sort = typeof sort !== 'undefined' ? sort : false;	
	genome = typeof genome !== 'undefined' ? genome : 'hg19';

	// generate json
	var json = generate_json( data, filter, sort );

	// if too many tracks, warning
	if ( json.length < 4 ) {
		alert("No data selected to visualize");
		return;
	}
	else if ( json.length > track_num_alert )
		if ( !confirm("The genome browser will visualize more than " + track_num_alert + " tracks ("+json.length+").\nAre you sure to proceed?" ) ) return;

	// save json in file system and send it to genome browser
	var req= new XMLHttpRequest();

	req.onreadystatechange = function() {
		if(req.readyState==4 && req.status==200) {
			var t=req.responseText;
			if(t.substr(0,5)=='ERROR') {
				alert('Failed to post data to server (error in cgi-bin app?)'+t,3);
			} else {
				var url = wubrowser_url + "?genome="+genome+"&tknamewidth="+width_trackname+"&datahub="+tmp_json_url+t;

				_DEBUG_&&console.log(tmp_json_url+t);

				if ( embed_id == null ) {
					window.open( url, "_blank");
				}
				else {					
					embed_wubrowser( embed_id, tmp_json_url+t, genome );
				}
			}
		}
	};

	var json_str = JSON.stringify(json);

	// save json at tmp_json_url
	req.open('POST', 'https://egg2.wustl.edu/cgi-bin/writejson', true);
	req.setRequestHeader("Content-type","application/x-www-form-urlencoded");
	req.send(json_str);
}

function get_chr_state_categories( chr_state ) 
{
	var obj_categories = { "-1": ["no information", "transparent"] }; 

	for(var i = 0; i < chr_state.length; i++ ) {
		var st = chr_state[i];
		var stateno = st.stateno;

		obj_categories[ stateno.toString() ] = [ st.desc, rgb_to_hex( st.rgb ) ];
	}

	return obj_categories;
}


function generate_json( data, filter, sort ) // data is an array of json objects
{
	sort = typeof sort !== 'undefined' ? sort : false;	

	// read default json
	var json = JSON.parse(JSON.stringify(data_wubrowser_json_base));

	// temp var for faster search
	var epg_json = new Object;
	var assay_json = new Object;
	
	// facet object
	var facet = { 
		"facet_table" :	[
			"Sample",
			"Assay",
			"Subtype",
		],
		"show_terms" : {
			"Sample": [
				"Sample",
			],
			"Assay": [
				"Assay",
			],
			"Subtype": [
				"Subtype",
			],
		},
		"type" : "metadata",
		"vocabulary_set" : { 
			"Sample" : {
				"terms" : {
				},
				"vocabulary" : {
					"Sample": {
						"Epigenomes" : [],
					},
				},
			},
			"Assay" : {
				"terms" : {
				},
				"vocabulary" : {
					"Assay": {
						"Assays" : [],
					},
				},
			},
			"Subtype" : {
				"terms" : {
				},
				"vocabulary" : {
					"Subtype": {
						"Categories" : [],
					},
				},
			},
		}
	};

	// add sample (epigenome) to facet table
	for(var k = 0; k < data_epg.length; k++ ) {
		var elem = data_epg[k];
		var order_str = elem.order.toString();
		var voc = facet.vocabulary_set.Sample;

		voc.terms[ order_str ] =  [ elem.mnemonic, ""];//, elem.color ];
		voc.vocabulary.Sample.Epigenomes.push( order_str );

		var order_str2 = (elem.groupid+20000).toString();
		var voc2 = facet.vocabulary_set[ "Subtype" ];

		voc2.terms[ order_str2 ] =  [ elem.group, "", elem.color ];
		voc2.vocabulary[ "Subtype" ].Categories.push( order_str2 );
	}

	// add assay to facet table
	for(var k = 0; k < data_assay.length; k++ ) {
		var elem = data_assay[k];
		var order_str = (elem.order+10000).toString(); // add 10000 to make it unique
		var voc = facet.vocabulary_set.Assay;

		voc.terms[ order_str ] =  [ elem.name, "", rgb_to_hex( elem.rgb ) ];//, elem.color ];
		voc.vocabulary.Assay.Assays.push( order_str );
	}
	
	json.push( facet );

	// temp array for faster search
	for(var k = 0; k < data_epg.length; k++ ) {
		var elem = data_epg[k];
		epg_json[ elem.eid ] = elem;		
	}

	for(var k = 0; k < data_assay.length; k++ ) {
		var elem = data_assay[k];
		assay_json[ elem.name ] = elem;
	}

	var has_align_data = false;

	// json object for each data set

	for(var k = 0; k < data.length; k++ ) {

		var json_current = new Array(); // store data with eid, will be sorted later

		for(var i = 0; i < data[k].data.length; i++ ) {
			data_elem = data[k].data[i];

			var obj = null;

			switch ( data[k].type ) {

				case "bigwig": 
					// if track height is specified in json obj
					var height_data = data[k]["height"];
					var height_row_bigwig = height_data==null ? height_row_bigwig_default : height_data;

					// parse rgb
					var arr = null;
					if ( data_elem.assay )	{
						if ( assay_json[ data_elem.assay ].rgb ) {
							arr = assay_json[ data_elem.assay ].rgb.split(",");
						}
						else {
							arr = [255,0,0];
						}
					}
					else {
						arr = data[k].rgb.split(","); 
					}

					var pr = parseInt( arr[0] ); var pg = parseInt( arr[1] ); var pb = parseInt( arr[2] );

					obj = { "type": "bigwig", "name": "", "url": "", "mode": 1,
						"qtc" : { "height": height_row_bigwig, "summeth":2, 
							"pr": pr, "pg": pg, "pb": pb, "smooth":3 } }; // pr pg pb: rgb color

					obj.qtc[ "thtype" ]= 1;


					if ( data_elem.assay && assay_json[ data_elem.assay ].name == "DNAMethylSBS" )	{
						obj.qtc[ "thmin" ] = 0;
						obj.qtc[ "thmax" ] = 1;

						if ( data[k].hasOwnProperty( "thtype" ) )
							obj.qtc[ "thtype" ] = data[k].thtype;
						if ( data[k].hasOwnProperty( "thmin" ) )
							obj.qtc[ "thmin" ] = data[k].thmin;
						if ( data[k].hasOwnProperty( "thmax" ) )
							obj.qtc[ "thmax" ] = data[k].thmax;
					}
					else if ( data_elem.assay && assay_json[ data_elem.assay ].signal == "L" )	{
						obj.qtc[ "thmin" ] = 2;
						obj.qtc[ "thmax" ] = 15;
					}
					else {
						obj.qtc[ "thmin" ] = 2;
						obj.qtc[ "thmax" ] = 40;
					}			
 
					if ( data_elem.info == "positive" ) {
						obj.qtc[ "thmin" ] = 0;
						obj.qtc[ "thmax" ] = 40;
					}
					else if ( data_elem.info == "negative" ) {
						obj.qtc[ "thmin" ] = -40;
						obj.qtc[ "thmax" ] = 0;
					}
			

					break;

				case "quantitativeCategorySeries":
	
					var height_row_qcat = height_row_qcat_default;

					obj= { "categories": [], "type": "quantitativeCategorySeries", 
						"name": "", "url": "", "mode": "show", "backgroundcolor":"#000000",
						"qtc": { "height" : height_row_qcat } };

					obj.qtc.height = data_elem.height;

					var chr_state = data[k].chr_state

					if ( chr_state ) {
						switch ( chr_state ) { 
							case "15":
								obj.categories = get_chr_state_categories( data_chr_state_15 );
								break;
							case "18":
								obj.categories = get_chr_state_categories( data_chr_state_18 );
								break;
							case "25":
								obj.categories = get_chr_state_categories( data_chr_state_25 ); 
								break;
							default : // if chr_state is embedded in data
								obj.categories = get_chr_state_categories( chr_state ); 
								break;
						}
					}
				
					break;

				case "categoryMatrix":

					var height_row_catmat = height_row_catmat_default;

					obj= { "categories": [], "type": "categoryMatrix", "name": "", "url": "", 
						"mode": "show", "rowcount":0, "rowheight":2, };
	
					obj.rowcount = data_elem.rowcount;
					obj.rowheight = data_elem.rowheight;

					var chr_state = data[k].chr_state

					if ( chr_state ) {
						switch ( chr_state ) { 
							case "15":
								obj.categories = get_chr_state_categories( data_chr_state_15 );
								break;
							case "18":
								obj.categories = get_chr_state_categories( data_chr_state_18 );
								break;
							case "25":
								obj.categories = get_chr_state_categories( data_chr_state_25 ); 
								break;
							default : // if chr_state is embedded in data
								obj.categories = get_chr_state_categories( chr_state ); 
								break;
						}
					}

					break;

				case "categorical":
					var height_row_cat = height_row_cat_default;


					obj= { "categories": [], "type": "categorical", "name": "", "url": "", "mode": 1,
						"qtc": { "height" : height_row_cat } };

					var chr_state = data[k][ "chr_state" ]

					if ( chr_state ) {
						switch ( chr_state ) { 
							case "15":
								obj.categories = get_chr_state_categories( data_chr_state_15 );
								break;
							case "18":
								obj.categories = get_chr_state_categories( data_chr_state_18 );
								break;
							case "25":
								obj.categories = get_chr_state_categories( data_chr_state_25 ); 
								break;
							default : // if chr_state is embedded in data
								obj.categories = get_chr_state_categories( chr_state ); 
								break;
						}
					}
				
					break;

				case "bed":
					var height_row_bed = height_row_bed_default;

					obj= { "type": "bed", "name": "", "url": "", "mode": 2, 
						"qtc": { "height" : height_row_bed, 
							"fontsize":"0pt","fontfamily":"sans-serif","fontbold":false, } };
				
					if ( data_elem.url_tail.search("tagAlign") >= 0 )
						has_align_data = true; // compact mode for tagAlign bed
					else if ( data_elem.url_tail.search("ReadCoverage.bed.gz") >= 0 )
						has_align_data = true;
					else if ( data_elem.url_tail.search("FractionalMethylation.bed.gz") >= 0 )
						has_align_data = true;
					else if ( data_elem.url_tail.search("Peak") >= 0 ) // imputed peak calls
						obj.mode = 6;
					else
						obj.mode = 1; // full mode for other beds

					break;

				case "interaction":
					var height_row_interact = height_row_interact_default;

					obj= { "type": "interaction", "name": "", "url": "", "mode": 4, 
						"qtc": { "anglescale":1,
							"fontsize":"10pt","fontfamily":"sans-serif","fontbold":false, } };
					break;

				case "hammock": 
					var height_row_hammock = height_row_hammock_default

					var arr = null;
					if ( data_elem.assay )	{
						if ( assay_json[ data_elem.assay ].rgb ) {
							arr = assay_json[ data_elem.assay ].rgb.split(",");
						}
						else {
							arr = [255,0,0];
						}
					}
					else {
			 			arr = data[k].rgb.split(","); 
					}

					var pr = parseInt( arr[0] ); var pg = parseInt( arr[1] ); var pb = parseInt( arr[2] );

					obj = { "type": "hammock", "name": "", "url": "", 
						"mode": "barplot", "showscoreidx":1, "boxcolor":rgb_to_hex( pr+","+pg+","+pb ),
						strokecolor:"#ff6600", 
						scorescalelst:[{"type":1,"min":0,"max":40},{"type":0},{"type":0},{"type":0}],
						scorenamelst:["signal value", "P value (-log10)","Q value (-log10)"],
							"qtc": { "height" : height_row_hammock, "summeth":2,
								"fontsize":"0pt","fontfamily":"sans-serif",
								 "nr": pr, "ng": pg, "nb": pb, } };

					// if plot type is defined (barplot, full, thin, dense)
					if ( data[k].mode != null ) {
						obj[ "mode" ] = "barplot"; //data[k].mode;
					}

					break;


				default :
					break;
			}
	
			obj.name = "";
			// for sorting by epg order
			if ( data_elem[ "eid" ] != null && data_elem[ "eid" ] != "" )
			        obj[ "epg_order" ] = epg_json[ data_elem.eid ].order;
			else
				obj[ "epg_order" ] = 10000;

			// for sorting by assay order
			if ( data_elem[ "assay" ] != null && data_elem[ "assay" ] != "" ) {
				obj[ "assay_order" ] = assay_json[ data_elem.assay ].order;
			}
			else
			        obj[ "assay_order" ] = 0;

			// for sorting RNA seq data (positive and then negative)
			if ( data[k].type == "interaction" ) {
				obj[ "etc_order" ] = 10000;
			}
			else if ( data_elem[ "info" ] != null ) {
				if ( data_elem[ "info" ] == "positive" )
					obj[ "etc_order" ] = 5;
				else if ( data_elem[ "info" ] == "negative" ) 
					obj[ "etc_order" ] = 6;
				else 
					obj[ "etc_order" ] = 0;
			}			
			else {
					obj[ "etc_order" ] = 0;
			}


			var epg_info = data_elem.eid ? 
				(data_elem.eid.length < 5 ? data_elem.eid + " " + epg_json[ data_elem.eid ].name + " " 
								: data_elem.eid + " " ) : "";

			var assay_info = data_elem.assay ? 
				data_elem.assay + " ": "";

			if ( data[k].type == "quantitativeCategorySeries" ) { 
				obj.name = data_elem.name;
			}
			else if ( data[k].type == "categoryMatrix" ) { 
				obj.name = data_elem.name;
			}
			else {
				for (var j=0; j< filter.length;j++ ) {
					var filt_elem = filter[j];

					if ( filt_elem.ignore_uc && data[k].uc )
						break;

					if ( data[k].type == "interaction" ) {
						if ( data_elem.eid == "" )
							obj.name = "All EIDs ";
						if ( filt_elem.eid == data_elem.eid )
							obj.name = epg_info;
					}
					else if ( data[k].type == "categorical" && filt_elem.eid == data_elem.eid ) { 
						// for chromatin state map, check eid only

						if ( filt_elem.assay == assay_name_chr_prim15  && data[k].chr_state == "15" )
							obj.name = epg_info;

						if ( filt_elem.assay == assay_name_chr_exp18   && data[k].chr_state == "18" )
							obj.name = epg_info;

						if ( filt_elem.assay == assay_name_imp_chr_12m && data[k].chr_state == "25" )
							obj.name = epg_info;

						if ( !filt_elem.assay && data[k].chr_state == "15" ) 
							obj.name = epg_info;
							//console.log(obj.name + " : " + data[k].chr_state + " : " + epg_info)
					}
					else if ( !filt_elem.eid && !filt_elem.assay ) { // no filter show all data
						obj.name = assay_info + epg_info;
					}					
					else if ( !filt_elem.eid && filt_elem.assay == data_elem.assay ) {
						obj.name = epg_info;
					}
					else if ( !filt_elem.assay && filt_elem.eid == data_elem.eid ) {
						obj.name = assay_info;
					}
					else if ( filt_elem.eid == data_elem.eid && filt_elem.assay == data_elem.assay ) {
						obj.name = assay_info + epg_info;
					}
					
					if ( data_elem.info.search("cluster") >= 0 ) {
						obj.name = " "
					}	

					if ( data[k].type == "hammock" && filt_elem.dpeak && obj.name != "" ) { 
						// if data-specific default peak should be used.

						//alert( filt_elem.dpeak  );
						var peak_type = assay_json[ data_elem.assay ].peak;
			
						//alert (peak_type );
						if ( data[k].info == "Narrow peaks" && ( peak_type == "N" || peak_type == "?" ) ) {
						}
						else if ( data[k].info == "Broad peaks" && ( peak_type == "B" ) ) {
						}
						else if ( data[k].info == "Gapped peaks" && ( peak_type == "G" ) ) {
						}
						else {
							obj.name = "";
						}
					}

					if ( obj.name != "" ) break;
				}
			}

			if ( obj.name != "" ) { // filter out if no obj name


				// add metadata
				var metadata = new Object;

				if ( data_elem.eid && data_elem.eid != "" ) { 
					metadata[ "Sample" ]= epg_json[ data_elem.eid ].order.toString();
					metadata[ "Subtype" ]= (epg_json[ data_elem.eid ].groupid + 20000).toString();
				}		
				if ( data_elem.assay && data_elem.assay != "" ) { 
					metadata[ "Assay" ]=  ( assay_json[ data_elem.assay ].order + 10000).toString(); 
					// +10000 to make ID unique
				}		

				if ( data_elem.info ) {
					obj.name = obj.name + "("+data_elem.info+") ";
				}

				if ( data.length > 1 && data[k].info.length ) { // if multiple data sets are used, specify data_info in track label
					obj.name = obj.name + "("+ data[k].info +") ";
				}

				// replace 'Input' with 'Control'
				obj.name = obj.name.replace("Input","Control");

				// replacece '_' with ' '
				obj.name = obj.name.replace(/_/g," ");

				obj[ "metadata" ] = metadata;
				obj.url  = data[k].url_head + data_elem.url_tail;

				json_current.push(obj);
			}
		}
		// sort data by epg order, assay order and info

		// Chrome array.sort is not stable so use the following.
		for(var i=0; i < json_current.length; i++) json_current[i]["position"] = i;
		json_current.sort( function (a, b) { if (a.etc_order === b.etc_order) 	return a.position - b.position; 	return a.etc_order-b.etc_order; } );

		for(var i=0; i < json_current.length; i++) json_current[i]["position"] = i;
		json_current.sort( function (a, b) { if (a.epg_order === b.epg_order) 	return a.position - b.position; 	return a.epg_order-b.epg_order; } );

		for(var i=0; i < json_current.length; i++) json_current[i]["position"] = i;
		json_current.sort( function (a, b) { if (a.assay_order === b.assay_order) return a.position - b.position;	return a.assay_order-b.assay_order; } );

		json = json.concat( json_current );
	}

	// special ordering of tracks for enh-gene interaction plot
	if ( sort ) {
		for(var i=0; i < json.length; i++) json[i]["position"] = i;
		json.sort( function (a, b) { if (a.etc_order === b.etc_order) 	return a.position - b.position; 	return a.etc_order-b.etc_order; } );

		for(var i=0; i < json.length; i++) json[i]["position"] = i;
		json.sort( function (a, b) { if (a.epg_order === b.epg_order) 	return a.position - b.position; 	return a.epg_order-b.epg_order; } );
	}

	// if there is alignment data, narrow down view
	if ( has_align_data ) {
		for(var i=0; i < json.length; i++) {
			if ( json[i].type && json[i].type == "coordinate_override" ) {
				//json[i].coord = "chr9,36329955,chr9,37537411" // original
				//json[i].coord = "chr9,36782751,chr9,37084615" // narrow
				json[i].coord = "chr9,36782751,chr9,36884615" // narrower
				break;
			}
		}
	}

	return json;
}

function populate_dropdown_by_eid( elem_id, data ) {

	var select_eid = document.getElementById( elem_id );
	var data_cnt = new Object;

	// loop through data, add eid to list
	for(var k = 0; k < data.length; k++ ) {
		var elem1 = data[k];

		for(var j=0; j < elem1.data.length; j++) {
			var elem2 = elem1.data[j];

			if ( data_cnt[ elem2.eid ] ) {
				data_cnt[ elem2.eid ]++;
			}
			else {
				data_cnt[ elem2.eid ] = 1;
			}	
		}
	}
	
	// populate
	for(var i = 0; i < data_epg.length; i++) {
		var elem = data_epg[i];

		if ( data_cnt[ elem.eid ] ) {
			var el = document.createElement("option");
			el.value = elem.eid;
			el.textContent = elem.eid.replace(/_/g," ") + " " + elem.name + " (" + data_cnt[ elem.eid ] + " trk)";
			select_eid.appendChild(el);
		}
	}	
}

function populate_dropdown_by_assay( elem_id, data ) {

	var select_assay = document.getElementById( elem_id );
	var data_cnt = new Object;

	// loop through data, add assay to list
	for(var k = 0; k < data.length; k++ ) {
		var elem1 = data[k];

		for(var j=0; j < elem1.data.length; j++) {
			var elem2 = elem1.data[j];

			if ( data_cnt[ elem2.assay ] ) {
				data_cnt[ elem2.assay ]++;
			}
			else {
				data_cnt[ elem2.assay ] = 1;
			}	
		}
	}

	// populate
	for(var i = 0; i < data_assay.length; i++) {
		var elem = data_assay[i];

		if ( data_cnt[ elem.name ] ) {
			
			var el = document.createElement("option");
			el.value = elem.name;
			el.textContent = elem.desc + " (" + data_cnt[ elem.name ] + " trk)";
			select_assay.appendChild(el);
		}
	}
}

function onchange_dropdown_by_eid(dropdown, data, paired_ddl_id, chkbox_id, embed_id, genome )
{
	if ( paired_ddl_id != "" ) {
		var select_assay = document.getElementById( paired_ddl_id );
		select_assay.selectedIndex = 0;
	}
	genome = typeof genome !== 'undefined' ? genome : 'hg19'

	var myindex  = dropdown.selectedIndex;
	var eid = dropdown.options[myindex].value;

	var embed_id = document.getElementById(chkbox_id).checked ? null : embed_id; // embed or null

	visualize_data( data, [{"eid":eid}], embed_id, false, genome );

    return true;
}

function onchange_dropdown_by_assay(dropdown, data, paired_ddl_id, chkbox_id, embed_id, genome )
{
	if ( paired_ddl_id != "" ) {
		var select_eid = document.getElementById( paired_ddl_id );
		select_eid.selectedIndex = 0;
	}
	genome = typeof genome !== 'undefined' ? genome : 'hg19'

	var myindex  = dropdown.selectedIndex;
	var assay = dropdown.options[myindex].value;

	var embed_id = document.getElementById( chkbox_id ).checked ? null : embed_id; // embed or null

	visualize_data( data, [{"assay":assay}], embed_id, false, genome );

	return true;
}

function onclick_btn_rna_bigwig( btn )
{
	var embed_id = document.getElementById("chkbox_new_page_rna_bigwig").checked ? null : "embed_rna_bigwig"; // embed or null
	visualize_data( [data_rna_bigwig], [{"":""}], embed_id );
	return true;
}

function make_random_str( len )
{
    var text = "";
    var possible = "0123456789";

    for( var i=0; i < len; i++ )
        text += possible.charAt(Math.floor(Math.random() * possible.length));

    return text;
}

function sleep( ms )
{
	setTimeout(function(){}, ms );
}

function init_rep_hist()
{
	// populate dropdown list

	// consolidated epg
	populate_dropdown_by_eid(  "ddl_eid_rep_hist_bw", [data_rep_hist_bw]);
	populate_dropdown_by_assay("ddl_assay_rep_hist_bw",[data_rep_hist_bw]);
	populate_dropdown_by_eid(  "ddl_eid_rep_hist_bw2", [data_rep_hist_bw2]);
	/*populate_dropdown_by_eid(  "ddl_eid_rep_hist_npeak"	, [data_rep_hist_npeak]);
	populate_dropdown_by_assay("ddl_assay_rep_hist_npeak"	, [data_rep_hist_npeak]);
	populate_dropdown_by_eid(  "ddl_eid_rep_hist_bpeak"	, [data_rep_hist_bpeak]);
	populate_dropdown_by_assay("ddl_assay_rep_hist_bpeak"	, [data_rep_hist_bpeak]);*/
}

function onchange_ddl_assay_rep_hist_bw(dropdown)
{
	return onchange_dropdown_by_assay(dropdown, [data_rep_hist_bw], "ddl_eid_rep_hist_bw", "chkbox_new_page_bw", "embed_rep_hist_bw" );
}

function onchange_ddl_eid_rep_hist_bw(dropdown)
{
	return onchange_dropdown_by_eid(dropdown, [data_rep_hist_bw], "ddl_assay_rep_hist_bw", "chkbox_new_page_bw", "embed_rep_hist_bw" );
}

function onchange_ddl_eid_rep_hist_bw2(dropdown)
{
	return onchange_dropdown_by_eid(dropdown, [data_rep_hist_bw2], "ddl_assay_rep_hist_bw2", "chkbox_new_page_bw", "embed_rep_hist_bw" );
}

function onchange_ddl_assay_rep_hist_npeak(dropdown)
{
	return onchange_dropdown_by_assay(dropdown, [data_rep_hist_npeak], "ddl_eid_rep_hist_npeak", "chkbox_new_page_npeak", "embed_rep_hist_npeak" );
}

function onchange_ddl_eid_rep_hist_npeak(dropdown)
{
	return onchange_dropdown_by_eid(dropdown, [data_rep_hist_npeak], "ddl_assay_rep_hist_npeak", "chkbox_new_page_npeak", "embed_rep_hist_npeak" );
}

function onchange_ddl_assay_rep_hist_bpeak(dropdown)
{
	return onchange_dropdown_by_assay(dropdown, [data_rep_hist_bpeak], "ddl_eid_rep_hist_bpeak", "chkbox_new_page_bpeak", "embed_rep_hist_bpeak" );
}

function onchange_ddl_eid_rep_hist_bpeak(dropdown)
{
	return onchange_dropdown_by_eid(dropdown, [data_rep_hist_bpeak], "ddl_assay_rep_hist_bpeak", "chkbox_new_page_bpeak", "embed_rep_hist_bpeak" );
}

function onchange_ddl_assay_rep_hist_gpeak(dropdown)
{
	return onchange_dropdown_by_assay(dropdown, [data_rep_hist_gpeak], "ddl_eid_rep_hist_gpeak", "chkbox_new_page_gpeak", "embed_rep_hist_gpeak" );
}

function onchange_ddl_eid_rep_hist_gpeak(dropdown)
{
	return onchange_dropdown_by_eid(dropdown, [data_rep_hist_gpeak], "ddl_assay_rep_hist_gpeak", "chkbox_new_page_gpeak", "embed_rep_hist_gpeak" );
}

/////////////////////////////// index.html (main page, grid visualization)

function init_vis()
{
	// checkboxes
	document.getElementById("chkbox_vis_rep_hist_pval").checked 	= true;
	document.getElementById("chkbox_vis_rep_hist_fc").checked 	= false;
	document.getElementById("chkbox_vis_imp_sig").checked 		= false;

	document.getElementById("chkbox_vis_rep_hist_dpeak").checked 	= true; // default
	document.getElementById("chkbox_vis_rep_hist_dpeak").disabled 	= false; //removeAttribute( "disabled", "false" );

		document.getElementById("chkbox_vis_rep_hist_npeak").checked 	= false;
		document.getElementById("chkbox_vis_rep_hist_npeak").disabled 	= true; 

		//document.getElementById("chkbox_vis_rep_hist_bpeak").checked 	= false;
		//document.getElementById("chkbox_vis_rep_hist_bpeak").disabled 	= true; 

		document.getElementById("chkbox_vis_rep_hist_gpeak").checked 	= true;
		document.getElementById("chkbox_vis_rep_hist_gpeak").disabled 	= true; 

	document.getElementById("chkbox_vis_rep_hist_align").checked 	= false;
	document.getElementById("chkbox_vis_chr_state_15_catMat").checked = false;
/*
	document.getElementById("chkbox_vis_chr_state_18_catMat").checked = false;
	document.getElementById("chkbox_vis_chr_state_25_catMat").checked = false;
*/
	document.getElementById("chkbox_vis_DMRs").checked = false;
	document.getElementById("chkbox_vis_RRBS").checked = false;
	document.getElementById("chkbox_vis_WGBS").checked = false;
	document.getElementById("chkbox_vis_mCRF").checked = false;

	// track height
	document.getElementById("text_vis_rep_hist_pval").value = height_row_bigwig_default;
	document.getElementById("text_vis_rep_hist_fc").value 	= height_row_bigwig_default;
	document.getElementById("text_vis_imp_sig").value 	= height_row_bigwig_default;
}

function visualize_grid( filter ) //filter = array of [{"eid":eid, "assay":assay }]
{
	// check if track height input is valid
	var height_vis_rep_hist_pval = parseInt( document.getElementById("text_vis_rep_hist_pval").value );
	var height_vis_rep_hist_fc   = parseInt( document.getElementById("text_vis_rep_hist_fc").value );
	var height_vis_imp_sig 	     = parseInt( document.getElementById("text_vis_imp_sig").value );

	if ( isNaN( height_vis_rep_hist_pval ) ) {
		alert("Invalid track height for uniform signal tracks (p-val)");
		return;
	}
	if ( isNaN( height_vis_rep_hist_fc ) ) {
		alert("Invalid track height for uniform signal tracks (fold enrichment)");
		return;
	}
	if ( isNaN( height_vis_imp_sig ) ) {
		alert("Invalid track height for imputed signal tracks");
		return;
	}

	if ( 	height_vis_rep_hist_pval < 10 || height_vis_rep_hist_pval > 40 ||
		height_vis_rep_hist_fc < 10 || height_vis_rep_hist_fc > 40 ||
		height_vis_imp_sig < 10 || height_vis_imp_sig > 40 ) {
		alert("Input track height between 10 and 40");
		return;
	}

	if ( document.getElementById("chkbox_vis_rep_hist_dpeak").checked ) {
		filter[ "peak" ] = true;
	}

	// set track height
	data_rep_hist_pval[ "height" ] 	= height_vis_rep_hist_pval;
	data_rep_hist_fc[ "height" ] 	= height_vis_rep_hist_fc;
	data_imp_sig[ "height" ] 	= height_vis_imp_sig;

	var data = concat_data();

	var embed_id = document.getElementById("chkbox_vis_new_page").checked ? null : "embed_vis"; // embed or null

	visualize_data( data, filter, embed_id );
}

function concat_data()	// prepare to gather data to visualize
{
	var data = new Array();

	data.push( data_imp_chr_12m );
	data.push( data_chr_prim15 );
	data.push( data_chr_exp18 );
/*	
	if ( document.getElementById("chkbox_vis_chr_state_25_catMat").checked ) {
		data.push( data_imp_chr_12m_catMat );
	}
*/
	if ( document.getElementById("chkbox_vis_chr_state_15_catMat").checked ) {
		data.push( data_chr_prim15_catMat );
	}
/*
	if ( document.getElementById("chkbox_vis_chr_state_18_catMat").checked ) {
		data.push( data_chr_exp18_catMat );
	}
*/
	if ( document.getElementById("chkbox_vis_chr_state_15_qcat").checked ) {
		data.push( data_chr_prim15_qcat );
	}

	if ( document.getElementById("chkbox_vis_imp_sig").checked ) {
		data.push( data_imp_sig );
	}
	if ( document.getElementById("chkbox_vis_rep_hist_align").checked ) {
		data.push( data_rep_hist_align );
		data.push( data_rep_hist_align_uc );
	}
	if ( document.getElementById("chkbox_vis_rep_hist_pval").checked ) {
		data.push( data_rep_hist_pval );
		data.push( data_rep_hist_pval_uc );
		data.push( data_rna_bigwig );
	}
	if ( document.getElementById("chkbox_vis_rep_hist_fc").checked ) {
		data.push( data_rep_hist_fc );
		data.push( data_rep_hist_fc_uc );
	}

	// if default peak is selected concat all peak data
	if ( document.getElementById("chkbox_vis_rep_hist_dpeak").checked ) {
		data.push( data_rep_hist_npeak );
		data.push( data_rep_hist_npeak_uc );
		//data.push( data_rep_hist_bpeak );
		data.push( data_rep_hist_gpeak );
		data.push( data_rep_hist_gpeak_uc );
	}
	else {
		if ( document.getElementById("chkbox_vis_rep_hist_npeak").checked ) {
			data.push( data_rep_hist_npeak );
			data.push( data_rep_hist_npeak_uc );
		}
/*
		if ( document.getElementById("chkbox_vis_rep_hist_bpeak").checked ) {
			data.push( data_rep_hist_bpeak );
			data.push( data_rep_hist_bpeak_uc );
		}
*/
		if ( document.getElementById("chkbox_vis_rep_hist_gpeak").checked ) {
			data.push( data_rep_hist_gpeak );
			data.push( data_rep_hist_gpeak_uc );
		}
	}

	if ( document.getElementById("chkbox_vis_DMRs").checked ) {
		data.push( data_dnamethyl_DMR_RRBS )
		data.push( data_dnamethyl_DMR_WGBS )
	}

	if ( document.getElementById("chkbox_vis_RRBS").checked ) {
		data.push( data_dnamethyl_RRBS_FM ) 
		data.push( data_dnamethyl_RRBS_RC ) 
	}

	if ( document.getElementById("chkbox_vis_WGBS").checked ) {
		data.push( data_dnamethyl_WGBS_FM ) 
		data.push( data_dnamethyl_WGBS_RC ) 
	}

	if ( document.getElementById("chkbox_vis_mCRF").checked ) {
		data.push( data_dnamethyl_mCRF_FM ) 
	}

	return data;
}

function get_data_avail_map()
{
	var data = concat_data();

	var epg_json = new Object; // temp array for fast search

	var avail = new Object; // availibility matrix (epg by assay)

	var dpeak_checked = document.getElementById("chkbox_vis_rep_hist_dpeak").checked;

	// temp array for fast epg,assay info search
	var assay_json = new Object;
	for(var i = 0; i < data_assay.length; i++ ) {
		var elem1 = data_assay[i];
		assay_json[ elem1.name ] = elem1;
	}

	for(var i = 0; i < data_epg.length; i++ ) {
		var elem1 = data_epg[i];
		epg_json[ elem1. eid ] = elem1;

		avail[ elem1.eid ] = new Object;
		for(var j = 0; j < data_assay.length; j++ ) {
			var elem2 = data_assay[j];
			avail[ elem1.eid ][ elem2.name ] = { "consolidated" : 0, "unconsolidated" : 0 };
		}
	}
	
	// iterate through all data sets
	for(var k = 0; k < data.length; k++ ) {
		var elem1 = data[k];

		for(var i = 0; i < elem1.data.length; i++ ) {
			elem2 = elem1.data[i];

			if ( elem2.eid && elem2.eid != "" && elem2.assay && elem2.assay != "" ) {

				var cnt = 0;
				// if default peak should be used, count number of data carefully according to assay type
				if ( dpeak_checked ) { 
					var peak_type = assay_json[ elem2.assay ].peak;

					if ( elem1.info == "Narrow peaks" ){					
						if ( peak_type == "N" || peak_type == "?" ) 
							cnt++;
					}
					else if ( elem1.info == "Broad peaks" ) {
						if ( peak_type == "B" ) 
							cnt++;
					}
					else if ( elem1.info == "Gapped peaks" ) {
						if ( peak_type == "G" ) 
							cnt++;
					}
					else {
						cnt++
					}
				}
				else { 
					cnt++;
				}

				if ( elem1.uc ) // if data is unconsolidated
					avail[ elem2.eid ][ elem2.assay ].unconsolidated += cnt;
				else
					avail[ elem2.eid ][ elem2.assay ].consolidated += cnt;
			}
		}
	}

	return avail;
}
