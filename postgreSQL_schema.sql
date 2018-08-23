--
-- PostgreSQL database dump
--

-- Dumped from database version 9.6.9
-- Dumped by pg_dump version 9.6.9

SET statement_timeout = 0;
SET lock_timeout = 0;
SET idle_in_transaction_session_timeout = 0;
SET client_encoding = 'UTF8';
SET standard_conforming_strings = on;
SELECT pg_catalog.set_config('search_path', '', false);
SET check_function_bodies = false;
SET client_min_messages = warning;
SET row_security = off;

--
-- Name: plpgsql; Type: EXTENSION; Schema: -; Owner: 
--

CREATE EXTENSION IF NOT EXISTS plpgsql WITH SCHEMA pg_catalog;


--
-- Name: EXTENSION plpgsql; Type: COMMENT; Schema: -; Owner: 
--

COMMENT ON EXTENSION plpgsql IS 'PL/pgSQL procedural language';


--
-- Name: make_uid(); Type: FUNCTION; Schema: public; Owner: username
--

CREATE FUNCTION public.make_uid() RETURNS text
    LANGUAGE plpgsql
    AS $$
BEGIN 
	return  CAST( '' || trunc(random()*10) ||  trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10)
                                || '-' || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10)
                                || '-' || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10)
                                || '-' || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10) || trunc(random()*10)
                                  AS varchar);
END;
$$;


ALTER FUNCTION public.make_uid() OWNER TO username;

SET default_tablespace = '';

SET default_with_oids = false;

--
-- Name: db_sequences; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.db_sequences (
    id text,
    upload_id integer,
    accession text,
    protein_name text,
    description text,
    sequence text,
    is_decoy boolean
);


ALTER TABLE public.db_sequences OWNER TO username;

--
-- Name: modifications; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.modifications (
    id bigint,
    upload_id integer,
    mod_name text,
    mass double precision,
    residues text,
    accession text
);


ALTER TABLE public.modifications OWNER TO username;

--
-- Name: peptide_evidences; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.peptide_evidences (
    upload_id integer,
    peptide_ref text,
    dbsequence_ref text,
    protein_accession text,
    pep_start integer,
    is_decoy boolean
);


ALTER TABLE public.peptide_evidences OWNER TO username;

--
-- Name: peptides; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.peptides (
    id text,
    upload_id integer,
    seq_mods text,
    link_site integer,
    crosslinker_modmass double precision,
    crosslinker_pair_id character varying
);


ALTER TABLE public.peptides OWNER TO username;

--
-- Name: protocols; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.protocols (
    id text NOT NULL,
    upload_id integer,
    protocol json,
    ms2_tol double precision
);


ALTER TABLE public.protocols OWNER TO username;

--
-- Name: spectra; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.spectra (
    id bigint,
    upload_id integer,
    peak_list text,
    peak_list_file_name text,
    scan_id text,
    frag_tol text,
    spectrum_ref text
);


ALTER TABLE public.spectra OWNER TO username;

--
-- Name: spectrum_identifications; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.spectrum_identifications (
    id bigint,
    upload_id integer,
    spectrum_id bigint,
    pep1_id text,
    pep2_id text,
    charge_state integer,
    pass_threshold boolean,
    rank integer,
    ions text,
    scores json,
    exp_mz double precision,
    calc_mz double precision,
    meta1 character varying,
    meta2 character varying,
    meta3 character varying
);


ALTER TABLE public.spectrum_identifications OWNER TO username;

--
-- Name: uploads; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.uploads (
    id integer NOT NULL,
    user_id integer,
    filename text,
    peak_list_file_names json,
    analysis_software json,
    provider json,
    audits json,
    samples json,
    analyses json,
    protocol json,
    bib json,
    spectra_formats json,
    upload_time date,
    default_pdb text,
    contains_crosslinks boolean,
    upload_error text,
    error_type text,
    upload_warnings json,
    origin text,
    random_id character varying DEFAULT public.make_uid()
);


ALTER TABLE public.uploads OWNER TO username;

--
-- Name: uploads_id_seq; Type: SEQUENCE; Schema: public; Owner: username
--

CREATE SEQUENCE public.uploads_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.uploads_id_seq OWNER TO username;

--
-- Name: uploads_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: username
--

ALTER SEQUENCE public.uploads_id_seq OWNED BY public.uploads.id;


--
-- Name: user_in_group; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.user_in_group (
    user_id integer,
    group_id integer
);


ALTER TABLE public.user_in_group OWNER TO username;

--
-- Name: users; Type: TABLE; Schema: public; Owner: username
--

CREATE TABLE public.users (
    user_name character varying,
    password character varying,
    email character varying,
    max_aas integer,
    max_spectra integer,
    gdpr_token character varying,
    id integer NOT NULL,
    gdpr_timestamp time without time zone,
    ptoken character varying,
    ptoken_timestamp character varying,
    hidden boolean
);


ALTER TABLE public.users OWNER TO username;

--
-- Name: users_id_seq; Type: SEQUENCE; Schema: public; Owner: username
--

CREATE SEQUENCE public.users_id_seq
    START WITH 1
    INCREMENT BY 1
    NO MINVALUE
    NO MAXVALUE
    CACHE 1;


ALTER TABLE public.users_id_seq OWNER TO username;

--
-- Name: users_id_seq; Type: SEQUENCE OWNED BY; Schema: public; Owner: username
--

ALTER SEQUENCE public.users_id_seq OWNED BY public.users.id;


--
-- Name: uploads id; Type: DEFAULT; Schema: public; Owner: username
--

ALTER TABLE ONLY public.uploads ALTER COLUMN id SET DEFAULT nextval('public.uploads_id_seq'::regclass);


--
-- Name: users id; Type: DEFAULT; Schema: public; Owner: username
--

ALTER TABLE ONLY public.users ALTER COLUMN id SET DEFAULT nextval('public.users_id_seq'::regclass);


--
-- Name: protocols protocols_pkey; Type: CONSTRAINT; Schema: public; Owner: username
--

ALTER TABLE ONLY public.protocols
    ADD CONSTRAINT protocols_pkey PRIMARY KEY (id);


--
-- Name: uploads uploads_pkey; Type: CONSTRAINT; Schema: public; Owner: username
--

ALTER TABLE ONLY public.uploads
    ADD CONSTRAINT uploads_pkey PRIMARY KEY (id);


--
-- PostgreSQL database dump complete
--

