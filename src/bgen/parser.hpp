#include "bgen.hpp"
#include "ImputeclassDetailed.h"

#include <jlibrary.h>
#include <zstream.h>

// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ProbSetter {
	typedef std::vector< std::vector< double > > Data ;
	ProbSetter( Data* result ):
		m_result( result ),
		m_sample_i(0)
	{}
		
	// Called once allowing us to set storage.
	void initialise( std::size_t number_of_samples, std::size_t number_of_alleles ) {
		m_result->clear() ;
		m_result->resize( number_of_samples ) ;
	}
	
	// If present with this signature, called once after initialise()
	// to set the minimum and maximum ploidy and numbers of probabilities among samples in the data.
	// This enables us to set up storage for the data ahead of time.
	void set_min_max_ploidy( uint32_t min_ploidy, uint32_t max_ploidy, uint32_t min_entries, uint32_t max_entries ) {
		for( std::size_t i = 0; i < m_result->size(); ++i ) {
			m_result->at( i ).reserve( max_entries ) ;
		}
	}
	
	// Called once per sample to determine whether we want data for this sample
	bool set_sample( std::size_t i ) {
		m_sample_i = i ;
		// Yes, here we want info for all samples.
		return true ;
	}
	
	// Called once per sample to set the number of probabilities that are present.
	void set_number_of_entries(
		std::size_t ploidy,
		std::size_t number_of_entries,
		genfile::OrderType order_type,
		genfile::ValueType value_type
	) {
		assert( value_type == genfile::eProbability ) ;
		m_result->at( m_sample_i ).resize( number_of_entries ) ;
		m_entry_i = 0 ;
	}

	// Called once for each genotype (or haplotype) probability per sample.
	void set_value( uint32_t, double value ) {
		m_result->at( m_sample_i ).at( m_entry_i++ ) = value ;
	}

	// Ditto, but called if data is missing for this sample.
	void set_value( uint32_t, genfile::MissingValue value ) {
		// Here we encode missing probabilities with -1
		m_result->at( m_sample_i ).at( m_entry_i++ ) = -1 ;
	}

	// If present with this signature, called once after all data has been set.
	void finalise() {
		// nothing to do in this implementation.
	}

private:
	Data* m_result ;
	std::size_t m_sample_i ;
	std::size_t m_entry_i ;
} ;


// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct PRSSetter {
    typedef std::vector<double> Data ;
    PRSSetter( Data* result , double factor, double rate, double mean):
        m_result( result ),
        m_factor(factor),
        m_rate(rate),
        m_mean(mean)
//        m_sample_i(0)
    {}

    // Called once allowing us to set storage.
    void initialise( std::size_t number_of_samples) {
        if(m_result->size()<number_of_samples){
            m_result->resize(number_of_samples);
        }
    }

    // Called once per sample to determine whether we want data for this sample
    bool set_sample( std::size_t i ) {
//        m_sample_i = i ;
        // Yes, here we want info for all samples.
        return true ;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_hap_value(int i, double value ) {
        double value2 = 1.0 - value;
        double meandiff = 0.5*(1-m_factor)+m_factor*value2-0.5*m_mean;
        (*m_result)[i] += meandiff * m_rate;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_genotype_value(int i, double value1, double value2 ) {
        double value = 2.0*max(1-value1-value2,0.0)+value2;
        double meandiff = 1-m_factor + m_factor*value - m_mean;
        (*m_result)[i] += meandiff * m_rate;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_ploidy1_genotype_value(int i, double value1 ) {
        double value = 1-m_factor + 2*m_factor*(1-value1) - m_mean;
        double meandiff = 1-m_factor + m_factor*value - m_mean;
        (*m_result)[i] += meandiff * m_rate;
    }

    void set_missing(int i){
        // do nothing
    }

    // If present with this signature, called once after all data has been set.
    void finalise() {
        // nothing to do in this implementation.
    }

private:
    Data* m_result ;
    double m_factor;
    double m_rate;
    double m_mean;
//    std::size_t m_sample_i ;
} ;


// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ImputeSetter {
    ImputeSetter( Imputeclass* result):
        m_result( result ),
        m_sqrsum(0.0),
        m_eafsum(0.0),
        m_psqrsum(0.0),
        m_ncount(0)
    {}

    // Called once allowing us to set storage.
    void initialise( std::size_t number_of_samples) {
        if(m_result->size()<number_of_samples){
            m_result->genotypes.resize(number_of_samples);
        }
        m_sqrsum=m_eafsum=m_psqrsum=0.0;
        m_ncount = 0;
    }

    // Called once per sample to determine whether we want data for this sample
    bool set_sample( std::size_t i ) {
//        m_sample_i = i ;
        // Yes, here we want info for all samples.
        return true ;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_hap_value(int i, double value ) {
        double w = 1.0 - value;
        (*m_result)[i] = w;
        m_sqrsum += w*w;
        m_eafsum += w;
        m_psqrsum += value;
        m_ncount += 1;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_genotype_value(int i, double value1, double value2 ) {
        double w = 2.0*max(1-value1-value2,0.0)+value2;
        (*m_result)[i] = w;
        m_sqrsum += w*w;
        m_eafsum += w;
        m_psqrsum += 4.0*max(1-value1-value2,0.0)+value2;
        m_ncount += 1;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_ploidy1_genotype_value(int i, double value1 ) {
        double w = 2*(1-value1);
        (*m_result)[i] = w;
        m_sqrsum += w*w;
        m_eafsum += w;
        m_psqrsum += 4*(1-value1);
        m_ncount += 1;
    }

    void set_missing(int i){
        (*m_result)[i] = -1.0;
    }

    // If present with this signature, called once after all data has been set.
    void finalise() {
        m_eafsum /= 2*m_ncount;
        m_result->eaf = m_eafsum;
        m_result->r2 = (m_sqrsum-4.0*m_ncount*m_eafsum*m_eafsum)/(m_psqrsum - 4.0*m_ncount*m_eafsum*m_eafsum);
        if(m_eafsum == 0 || m_eafsum==1) m_result->r2 = -1;
    }

private:
    Imputeclass* m_result ;
    double m_sqrsum;
    double m_eafsum;
    double m_psqrsum;
    int m_ncount;
} ;


// ProbSetter is a callback object appropriate for passing to bgen::read_genotype_data_block() or
// the synonymous method of genfile::bgen::View. See the comment in bgen.hpp above
// bgen::read_genotype_data_block(), or the bgen wiki for a description of the API.
// The purpose of this object is to store genotype probability values in the desired
// data structure (which here is a vector of vectors of doubles).
struct ImputeFilter {
    ImputeFilter( Imputeclass* result, const vector<int>* include):
        m_result( result ),
        m_include(include),
        m_sqrsum(0.0),
        m_eafsum(0.0),
        m_psqrsum(0.0),
        m_ncount(0)
    {}

    // Called once allowing us to set storage.
    void initialise( std::size_t number_of_samples) {
        m_result->genotypes.resize(0);
        m_sqrsum=m_eafsum=m_psqrsum=0.0;
        m_ncount = 0;
    }

    // Called once per sample to determine whether we want data for this sample
    bool set_sample( std::size_t i ) {
//        m_sample_i = i ;
        // Yes, here we want info for all samples.
        return true ;
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_hap_value(int i, double value ) {
        if((*m_include)[i]==1){
            double w = 1.0 - value;
            m_result->genotypes.push_back(w);
            m_sqrsum += w*w;
            m_eafsum += w;
            m_psqrsum += value;
            m_ncount += 1;
        }
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_genotype_value(int i, double value1, double value2 ) {
        if((*m_include)[i]==1){
            double w = 2.0*max(1-value1-value2,0.0)+value2;
            m_result->genotypes.push_back(w);
            m_sqrsum += w*w;
            m_eafsum += w;
            m_psqrsum += 4.0*max(1-value1-value2,0.0)+value2;
            m_ncount += 1;
        }
    }

    // Called once for each genotype (or haplotype) probability per sample.
    void update_ploidy1_genotype_value(int i, double value1 ) {
        if((*m_include)[i]==1){
            double w = 2*(1-value1);
            m_result->genotypes.push_back(w);
            m_sqrsum += w*w;
            m_eafsum += w;
            m_psqrsum += 4*(1-value1);
            m_ncount += 1;
        }
    }

    void set_missing(int i){
        if((*m_include)[i]==1){
            m_result->genotypes.push_back(-1.0);
        }
    }

    // If present with this signature, called once after all data has been set.
    void finalise() {
        m_eafsum /= 2*m_ncount;
        m_result->eaf = m_eafsum;
        m_result->r2 = (m_sqrsum-4.0*m_ncount*m_eafsum*m_eafsum)/(m_psqrsum - 4.0*m_ncount*m_eafsum*m_eafsum);
        if(m_eafsum == 0 || m_eafsum==1) m_result->r2 = -1;
    }

private:
    Imputeclass* m_result ;
    const vector<int>* m_include;
    double m_sqrsum;
    double m_eafsum;
    double m_psqrsum;
    int m_ncount;
} ;


// BgenParser is a thin wrapper around the core functions in genfile/bgen/bgen.hpp.
// This class tracks file state and handles passing the right callbacks.
struct BgenParser {
	
	BgenParser( std::string const& filename ):
		m_filename( filename ),
		m_state( e_NotOpen ),
		m_have_sample_ids( false )
	{
        is_bgenformat = (filename.size()>=4 && filename.substr(filename.size()-4)=="bgen");

        // Open the stream
        if(is_bgenformat)
        {
            m_stream.reset(
                        new std::ifstream( filename, std::ifstream::binary )
                        ) ;
            if( !*m_stream ) {
                throw std::invalid_argument( filename ) ;
            }
        }
        else
        {
            m_zstream.reset(new zifstream(filename));
        }

		m_state = e_Open ;

        if(is_bgenformat)
        {
            // Read the offset, header, and sample IDs if present.
            genfile::bgen::read_offset( *m_stream, &m_offset ) ;
            genfile::bgen::read_header_block( *m_stream, &m_context ) ;
            if( m_context.flags & genfile::bgen::e_SampleIdentifiers ) {
                genfile::bgen::read_sample_identifier_block(
                            *m_stream, m_context,
                            [this]( std::string id ) { m_sample_ids.push_back( id ) ; }
                ) ;
                m_have_sample_ids = true ;
            }

            // Jump to the first variant data block.
            m_stream->seekg( m_offset + 4 ) ;

            // We keep track of state (though it's not really needed for this implementation.)
            m_state = e_ReadyForVariant ;
        }
	}

	std::ostream& summarise( std::ostream& o ) const {
		o << "BgenParser: bgen file ("
			<< ( m_context.flags & genfile::bgen::e_Layout2 ? "v1.2 layout" : "v1.1 layout" )
			<< ", "
			<< ( m_context.flags & genfile::bgen::e_CompressedSNPBlocks ? "compressed" : "uncompressed" ) << ")"
			<< " with " 
			<< m_context.number_of_samples << " " << ( m_have_sample_ids ? "named" : "anonymous" ) << " samples and "
			<< m_context.number_of_variants << " variants.\n" ;
		return o ;
	}

	std::size_t number_of_samples() const {
		return m_context.number_of_samples ;
	}

    std::size_t number_of_variants() const {
        return m_context.number_of_variants ;
    }

	// Report the sample IDs in the file using the given setter object
	// (If there are no sample IDs in the file, we report a dummy identifier).
	template< typename Setter >
	void get_sample_ids( Setter setter ) {
		if( m_have_sample_ids ) {
			for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
				setter( m_sample_ids[i] ) ;
			}
		} else {
			for( std::size_t i = 0; i < m_context.number_of_samples; ++i ) {
				setter( "(unknown_sample_" + std::to_string( i+1 ) + ")" ) ;
			}
		}
	}

	// Attempt to read identifying information about a variant from the bgen file, returning
	// it in the given fields.
	// If this method returns true, data was successfully read, and it should be safe to call read_probs()
	// or ignore_probs().
	// If this method returns false, data was not successfully read indicating the end of the file.
	bool read_variant(
		std::string* chromosome,
		uint32_t* position,
		std::string* rsid,
		std::vector< std::string >* alleles
	) {
		assert( m_state == e_ReadyForVariant ) ;
		std::string SNPID ; // read but ignored in this toy implementation
		
		if(
			genfile::bgen::read_snp_identifying_data(
				*m_stream, m_context,
				&SNPID, rsid, chromosome, position,
				[&alleles]( std::size_t n ) { alleles->resize( n ) ; },
				[&alleles]( std::size_t i, std::string const& allele ) { alleles->at(i) = allele ; }
			)
		) {
			m_state = e_ReadyForProbs ;
			return true ;
		} else {
			return false ;
		}
	}
	
	// Read genotype probability data for the SNP just read using read_variant()
	// After calling this method it should be safe to call read_variant() to fetch
	// the next variant from the file.
	void read_probs( std::vector< std::vector< double > >* probs ) {
		assert( m_state == e_ReadyForProbs ) ;
		ProbSetter setter( probs ) ;
		genfile::bgen::read_and_parse_genotype_data_block< ProbSetter >(
            *m_stream,
			m_context,
			setter,
			&m_buffer1,
			&m_buffer2
		) ;
		m_state = e_ReadyForVariant ;
	}

	// Ignore genotype probability data for the SNP just read using read_variant()
	// After calling this method it should be safe to call read_variant()
	// to fetch the next variant from the file.
	void ignore_probs() {
        genfile::bgen::ignore_genotype_data_block( *m_stream, m_context ) ;
		m_state = e_ReadyForVariant ;
    }

    // Read genotype probability data for the SNP just read using read_variant()
    // After calling this method it should be safe to call read_variant() to fetch
    // the next variant from the file.
    void UpdatePRS( std::vector<double>& prs, double factor, double rate, double mean) {
        assert( m_state == e_ReadyForProbs ) ;
        PRSSetter setter( &prs, factor, rate, mean ) ;
        genfile::bgen::read_and_parse_genotype_data_block_prs<PRSSetter>(
            *m_stream,
            m_context,
            setter,
            &m_buffer1,
            &m_buffer2
        ) ;
        m_state = e_ReadyForVariant ;
    }

template<class T>
    bool ReadImputeclass(Imputeclass& imputeclass, const T& include)
    {
        uint32_t position;

        std::vector< std::string > alleles ;

        vector<vector<double>> probs;

        if(!is_bgenformat)
        {
            return ::ReadImputeclass(*m_zstream, imputeclass, include);
        }

        if(!read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);

        ImputeFilter setter(&imputeclass, &include);

        genfile::bgen::read_and_parse_genotype_data_block_prs<ImputeFilter>(
            *m_stream,
            m_context,
            setter,
            &m_buffer1,
            &m_buffer2
        ) ;
        m_state = e_ReadyForVariant ;

//        read_probs( &probs ) ;

//        int numinclude = count(include.begin(),include.end(),1);

//        imputeclass.genotypes.resize(numinclude);

//        int index = 0;

//        double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

//        int ncount = 0;

//        forc(i, probs)
//        {
//            if(include[i]==1)
//            {
//                if(probs[i].size()==3)
//                {
//                    imputeclass[index] = 2*probs[i][2]+probs[i][1];
//                }
//                else if(probs[i].size()==2)
//                {
//                    imputeclass[index] = 2*probs[i][1];
//                }
//                else
//                {
//                    imputeclass[index] = -3;
//                }

//                if(imputeclass[index]>-0.5)
//                {
//                    sqrsum += Sqr(imputeclass[index]);

//                    p1sum += (probs[i].size()==3?probs[i][1]:0);

//                    p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

//                    ++ncount;

//                    eafsum += imputeclass[index];
//                }

//                index+=1;
//            }
//        }

//        eafsum /= 2*ncount;

//        imputeclass.eaf = eafsum;

//        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

//        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

        return true;
    }

template<class T>
    bool ReadImputeclassCompare(Imputeclass& imputeclass, const T& include)
    {
        uint32_t position;

        std::vector< std::string > alleles ;

        vector<vector<double>> probs;

        if(!is_bgenformat)
        {
            return ::ReadImputeclass(*m_zstream, imputeclass, include);
        }

        if(!read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);

//        ImputeFilter setter(&imputeclass, &include);

//        genfile::bgen::read_and_parse_genotype_data_block_prs<ImputeFilter>(
//            *m_stream,
//            m_context,
//            setter,
//            &m_buffer1,
//            &m_buffer2
//        ) ;
//        m_state = e_ReadyForVariant ;

        read_probs( &probs ) ;

        int numinclude = count(include.begin(),include.end(),1);

        imputeclass.genotypes.resize(numinclude);

        int index = 0;

        double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0, psqrsum=0;

        int ncount = 0;

        forc(i, probs)
        {
            if(include[i]==1)
            {
                if(probs[i].size()==3)
                {
                    imputeclass[index] = 2*probs[i][2]+probs[i][1];
                }
                else if(probs[i].size()==2)
                {
                    imputeclass[index] = 2*probs[i][1];
                }
                else
                {
                    imputeclass[index] = -3;
                }

                if(imputeclass[index]>-0.5)
                {
                    sqrsum += Sqr(imputeclass[index]);

                    p1sum += (probs[i].size()==3?probs[i][1]:0);

                    p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

                    psqrsum += (probs[i].size()==3?4*probs[i][2]+probs[i][1]:4*probs[i][1]);

                    ++ncount;

                    eafsum += imputeclass[index];
                }

                index+=1;
            }
        }

        eafsum /= 2*ncount;

        imputeclass.eaf = eafsum;

//        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(psqrsum - 4.0*ncount*eafsum*eafsum);

        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

        return true;
    }

    bool ReadAllImputeclass(Imputeclass& imputeclass)
    {
        if(!is_bgenformat)
        {
            return bool(*m_zstream >> imputeclass);
        }

        uint32_t position;

        std::vector< std::string > alleles ;

        vector<vector<double>> probs;


        if(!read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);

        ImputeSetter setter(&imputeclass);

        genfile::bgen::read_and_parse_genotype_data_block_prs<ImputeSetter>(
            *m_stream,
            m_context,
            setter,
            &m_buffer1,
            &m_buffer2
        ) ;
        m_state = e_ReadyForVariant ;


//        read_probs( &probs ) ;

//        imputeclass.genotypes.resize(probs.size());

//        double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

//        forc(i, probs)
//        {
//            if(probs[i].size()==3)
//            {
//                imputeclass[i] = 2*probs[i][2]+probs[i][1];
//            }
//            else if(probs[i].size()==2)
//            {
//                imputeclass[i] = 2*probs[i][1];
//            }
//            else
//            {
//                imputeclass[i] = -1;
//            }

//            if(imputeclass[i]>-0.5)
//            {
//                sqrsum += Sqr(imputeclass[i]);

//                p1sum += (probs[i].size()==3?probs[i][1]:0);

//                p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

//                eafsum += imputeclass[i];
//            }
//        }

//        int ncount = probs.size();

//        eafsum /= 2*ncount;

//        imputeclass.eaf = eafsum;

//        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

//        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

        return true;
    }

    bool ReadAllImputeclassCompare(Imputeclass& imputeclass)
    {
        if(!is_bgenformat)
        {
            return bool(*m_zstream >> imputeclass);
        }

        uint32_t position;

        std::vector< std::string > alleles ;

        vector<vector<double>> probs;


        if(!read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);

//        ImputeSetter setter(&imputeclass);

//        genfile::bgen::read_and_parse_genotype_data_block_prs<ImputeSetter>(
//                *m_stream,
//                m_context,
//                setter,
//                &m_buffer1,
//                &m_buffer2
//                ) ;
//        m_state = e_ReadyForVariant ;


            read_probs( &probs ) ;

            imputeclass.genotypes.resize(probs.size());

            double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0, psqrsum =0;

            forc(i, probs)
            {
                if(probs[i].size()==3)
                {
                    imputeclass[i] = 2*probs[i][2]+probs[i][1];
                }
                else if(probs[i].size()==2)
                {
                    imputeclass[i] = 2*probs[i][1];
                }
                else
                {
                    imputeclass[i] = -1;
                }

                if(imputeclass[i]>-0.5)
                {
                    sqrsum += Sqr(imputeclass[i]);

                    p1sum += (probs[i].size()==3?probs[i][1]:0);

                    p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

                    psqrsum += (probs[i].size()==3?4.0*probs[i][2]+probs[i][1]:4.0*probs[i][1]);

                    eafsum += imputeclass[i];
                }
            }

            int ncount = probs.size();

            eafsum /= 2*ncount;

            imputeclass.eaf = eafsum;

//            imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

            imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(psqrsum - 4.0*ncount*eafsum*eafsum);

            if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

            return true;
    }

    void ReadAllProbs(Imputeclass& imputeclass){

        vector<vector<double>> probs;

        read_probs( &probs ) ;

        imputeclass.genotypes.resize(probs.size());

        double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

        forc(i, probs)
        {
            if(probs[i].size()==3)
            {
                imputeclass[i] = 2*probs[i][2]+probs[i][1];
            }
            else if(probs[i].size()==2)
            {
                imputeclass[i] = 2*probs[i][1];
            }
            else
            {
                imputeclass[i] = -1;
            }

            if(imputeclass[i]>-0.5)
            {
                sqrsum += Sqr(imputeclass[i]);

                p1sum += (probs[i].size()==3?probs[i][1]:0);

                p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

                eafsum += imputeclass[i];
            }
        }

        int ncount = probs.size();

        eafsum /= 2*ncount;

        imputeclass.eaf = eafsum;

        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;
    }

    bool ReadImputeclassDetailed(ImputeclassDetailed& imputeclass)
    {
        if(!is_bgenformat)
        {
            return bool(*m_zstream >> imputeclass);
        }

        uint32_t position;

        std::vector< std::string > alleles ;

        vector<vector<double>> probs;


        if(!read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);


        read_probs( &probs ) ;

        imputeclass.genotypes.resize(probs.size());

        imputeclass.probs.resize(probs.size());


        double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

        forc(i, probs)
        {
            if(probs[i].size()==3)
            {
                imputeclass[i] = 2*probs[i][2]+probs[i][1];

                forl(j,3) imputeclass.probs[i][j] = probs[i][j];
            }
            else if(probs[i].size()==2)
            {
                imputeclass[i] = 2*probs[i][1];

                tie(imputeclass.probs[i][0], imputeclass.probs[i][1],imputeclass.probs[i][2]) = make_tuple(probs[i][0],0,probs[i][1]);
            }
            else
            {
                imputeclass[i] = -3;

                imputeclass.probs[i][0] = imputeclass.probs[i][1]  = imputeclass.probs[i][2] = -1;
            }

            if(imputeclass[i]>-0.5)
            {
                sqrsum += Sqr(imputeclass[i]);

                p1sum += (probs[i].size()==3?probs[i][1]:0);

                p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

                eafsum += imputeclass[i];
            }
        }

        int ncount = probs.size();

        eafsum /= 2*ncount;

        imputeclass.eaf = eafsum;

        imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

        if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

        return true;
    }

        void JumpToPosition(size_t position){
            m_stream->seekg(position);
            m_state = e_ReadyForVariant ;
        }

    std::unique_ptr<zifstream > m_zstream ;

private:
    bool is_bgenformat;
	std::string const m_filename ;
	std::unique_ptr< std::istream > m_stream ;

	// bgen::Context object holds information from the header block,
	// including bgen flags
	genfile::bgen::Context m_context ;

	// offset byte from top of bgen file.
	uint32_t m_offset ;

	// We keep track of our state in the file.
	// Not strictly necessary for this implentation but makes it clear that
	// calls must be read_variant() followed by read_probs() (or ignore_probs())
	// repeatedly.
	enum State { e_NotOpen = 0, e_Open = 1, e_ReadyForVariant = 2, e_ReadyForProbs = 3, eComplete = 4 } ;
	State m_state ;

	// If the BGEN file contains samples ids, they will be read here.
	bool m_have_sample_ids ;
	std::vector< std::string > m_sample_ids ;
	
	// Buffers, these are used as working space by bgen implementation.
	std::vector< genfile::byte_t > m_buffer1, m_buffer2 ;
} ;


template<class T>
bool ReadImputeclass(BgenParser& bgenParser, Imputeclass& imputeclass, const T& include, bool read_position=true)
{
    uint32_t position;

    std::vector< std::string > alleles ;

    vector<vector<double>> probs;

    if(read_position)
    {
        if(!bgenParser.read_variant( &imputeclass.chromosome, &position, &imputeclass.rsnum, &alleles )) return false;

        tie(imputeclass.ref, imputeclass.alt) = {alleles[0],alleles[1]};

        imputeclass.position = int(position);
    }

    bgenParser.read_probs( &probs ) ;

    int numinclude = count(include.begin(),include.end(),1);

    imputeclass.genotypes.resize(numinclude);

    int index = 0;

    double p1sum=0, p2sum=0, sqrsum=0, eafsum = 0;

    int ncount = 0;

    forc(i, probs)
    {
        if(include[i]==1)
        {
            if(probs[i].size()==3)
            {
                imputeclass[index] = 2*probs[i][2]+probs[i][1];
            }
            else if(probs[i].size()==2)
            {
                imputeclass[index] = 2*probs[i][1];
            }
            else
            {
                imputeclass[index] = -3;
            }

            if(imputeclass[index]>-0.5)
            {
                sqrsum += Sqr(imputeclass[index]);

                p1sum += (probs[i].size()==3?probs[i][1]:0);

                p2sum += (probs[i].size()==3?probs[i][2]:probs[i][1]);

                ++ncount;

                eafsum += imputeclass[index];
            }

            index+=1;
        }
    }

    eafsum /= 2*ncount;

    imputeclass.eaf = eafsum;

    imputeclass.r2 = (sqrsum-4.0*ncount*eafsum*eafsum)/(p1sum + 4*p2sum - 4.0*ncount*eafsum*eafsum);

    if(eafsum == 0 || eafsum==1) imputeclass.r2 = -1;

    return true;
}
