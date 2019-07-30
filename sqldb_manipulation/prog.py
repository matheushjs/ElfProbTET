import sqlite3
import random
from faker import Faker

fake = Faker()
fake.seed(72)

conn = sqlite3.connect("database.db")

conn.execute("DROP TABLE IF EXISTS example;");

conn.execute("CREATE TABLE example(name CHAR(120), desc TEXT, age INTEGER)")

for i in range(200):
    conn.execute("INSERT INTO example(name, desc, age) VALUES ('{}', '{}', '{}')".format(fake.name(), fake.text(), random.randint(1, 100)))
    conn.commit()

res = conn.execute("SELECT * FROM example WHERE name LIKE 'M%'")

for r in res:
    print(r)
