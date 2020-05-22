from sqlalchemy import create_engine, Column, Sequence, Integer, String
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import sessionmaker
engine = create_engine("sqlite:///D:\\BMSIS_YSP\\test_rdkit\\network.db")
Session = sessionmaker(bind=engine)
Base = declarative_base()

def execute_query(sql, params):
    """
    General function for executing a query and expecting no data
    in return.
    """
    conn = engine.raw_connection()
    cur = conn.cursor()
    try:
        cur.execute(sql, params)
        cur.close()
        conn.commit()
    except Exception as e:
        print(f"Error executing query: {sql}... {params}... {e}.")
    finally:
        conn.close()


class Molecule(Base):
    """
    Save parent molecules (no reactants) and also
    child (product) molecules with reactant information,
    as well as which generation the product was formed.
    """
    __tablename__ = 'molecule'
    rowid = Column(Integer, Sequence('rowid'), primary_key=True)
    
    reactant_1 = Column(Integer) # self-referential foreign key
    reactant_2 = Column(Integer) # self-referential foreign key
    smiles_str = Column(String) # unique identifier for molecule to translate between formats
    generation_formed = Column(Integer) # generation (loop) number at which product was first formed

# will create the db file if not existing with the above specified table & constraints
Base.metadata.create_all(engine)