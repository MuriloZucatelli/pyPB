class Employee:
    def __init__(self, id, name):
        self.id = id
        self.name = name
        self.address = None

class Address:
    def __init__(self, street, city, state, zipcode, street2=''):
        self.street = street
        self.street2 = street2
        self.city = city
        self.state = state
        self.zipcode = zipcode

    def __str__(self):
        lines = [self.street]
        if self.street2:
            lines.append(self.street2)
        lines.append(f'{self.city}, {self.state} {self.zipcode}')
        return '\n'.join(lines)
    
manager = employees.Manager(1, 'Mary Poppins', 3000)
manager.address = contacts.Address(
    '121 Admin Rd', 
    'Concord', 
    'NH', 
    '03301'
)