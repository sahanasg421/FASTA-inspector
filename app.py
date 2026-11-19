import pandas as pd
import streamlit as st
from fasta_inspector import parse_fasta, analyse, summarise, EXAMPLE

st.title("FASTA Inspector")

text = st.text_area("Paste FASTA sequences", EXAMPLE, height=300)
min_length = st.number_input("Minimum length", min_value=1, value=20)

sequences, repeated_ids = parse_fasta(text)

if sequences:
    count, shortest, longest, mean_length, gc = summarise(sequences)
    st.write(f"{count} sequences. Shortest {shortest}, longest {longest}")
    st.write(f"Mean length {mean_length:.1f}, overall GC {gc:.1f}%")
    st.write("Repeated IDs:", repeated_ids)

    rows = analyse(sequences, min_length)
    st.dataframe(pd.DataFrame(rows, columns=["ID", "Length", "GC %", "Notes"]))
else:
    st.info("Paste some FASTA text above.")
